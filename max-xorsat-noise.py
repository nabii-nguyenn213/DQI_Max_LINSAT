import argparse
import importlib.util
import sys
import types
import time
import os
import csv
from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
import matplotlib.pyplot as plt 
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument(
    "--source",
    default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "dqi_collision_feedforward.py"
    )
)
parser.add_argument("--ell", type=int, default=4)
parser.add_argument("--trajectories", type=int, default=20)
parser.add_argument(
    "--error-rates",
    default="0,0.0001,0.0005,0.001"
)
parser.add_argument("--seed", type=int, default=20261111)
parser.add_argument(
    "--csv",
    default="",
    help="Optional noise-results CSV"
)
args = parser.parse_args()

@dataclass
class Op:
    typ: str
    wires: tuple
    controls: tuple = ()
    values: tuple = ()
    theta: float = 0.0
    prep: tuple = ()

ops = []
ctl_stack = []

def emit(
    typ,
    wires,
    controls=(),
    values=(),
    theta=0.0,
    prep=()
):
    inherited = sum(
        (list(ctl) for ctl, val in ctl_stack), []
    )

    inherited_vals = sum(
        (list(val) for ctl, val in ctl_stack), []
    )

    op = Op(
        typ,
        tuple(wires),
        tuple(inherited) + tuple(controls),
        tuple(inherited_vals) + tuple(values),
        float(theta),
        tuple(prep)
    )

    assert not set(op.wires) & set(op.controls), (
        op,
        "target cannot be a control"
    )

    ops.append(op)


def ctrl(f, control=(), control_values=None):
    ctl = (
        tuple(control)
        if isinstance(control, (tuple, list))
        else (control,)
    )

    vals = tuple(
        int(z)
        for z in (
            control_values
            if control_values is not None
            else [1] * len(ctl)
        )
    )

    def new(*a, **kw):
        ctl_stack.append((ctl, vals))

        try:
            return f(*a, **kw)
        finally:
            ctl_stack.pop()

    return new


def adjoint(f):
    def run(*a, **kw):
        old = len(ops)

        f(*a, **kw)

        block = ops[old:]
        del ops[old:]

        for op in reversed(block):
            ops.append(
                Op(
                    op.typ,
                    op.wires,
                    op.controls,
                    op.values,
                    -op.theta,
                    op.prep
                )
            )

    return run


def gate_func(t):
    def fn(theta=None, wires=None, **kwargs):

        if t == "CNOT":
            emit(
                "X",
                (wires[1],),
                (wires[0],),
                (1,)
            )

        elif t == "CRY":
            emit(
                "RY",
                (wires[1],),
                (wires[0],),
                (1,),
                theta
            )

        elif t == "RY":
            emit("RY", (wires,), theta=theta)

        elif t == "StatePrep":
            emit(
                "PREP",
                tuple(wires),
                prep=tuple(float(x) for x in theta)
            )

        else:
            emit(t, (wires,))

    return fn


# ============================================================
# 3. PENNYLANE RECORDING SHIM
# ============================================================

mock = types.ModuleType("pennylane")

for typ in [
    "PauliX",
    "PauliZ",
    "RY",
    "CNOT",
    "CRY",
    "Hadamard"
]:
    fn = gate_func(
        {
            "PauliX": "X",
            "PauliZ": "Z",
            "Hadamard": "H"
        }.get(typ, typ)
    )

    setattr(mock, typ, fn)


mock.StatePrep = lambda state, wires, pad_with=0: emit(
    "PREP",
    tuple(wires),
    prep=tuple(state)
)

mock.ctrl = ctrl
mock.adjoint = adjoint
mock.device = lambda *args, **kwargs: None
mock.qnode = lambda device: (lambda f: f)
mock.probs = lambda wires=None: None
mock.counts = lambda wires=None: None
mock.sample = lambda wires=None: None

sys.modules["pennylane"] = mock

argv_save = sys.argv

sys.argv = [
    args.source,
    "--ell",
    str(args.ell)
]

spec = importlib.util.spec_from_file_location(
    "dqi_imported",
    args.source
)

module = importlib.util.module_from_spec(spec)

spec.loader.exec_module(module)

sys.argv = argv_save


# ============================================================
# 4. CONSTRUCT LOGICAL CIRCUITS
# ============================================================

def make_ops(method, grover=0):
    ops.clear()

    module.dqi_circuit(method)

    for _ in range(grover):

        # Good-subspace phase oracle
        for w in module.error_wires:
            mock.PauliZ(wires=w)

        # Inverse state preparation
        mock.adjoint(module.dqi_circuit)("eraser")

        # Reflection about initial state
        module.reflect_initial_zero()

        # Forward state preparation
        module.dqi_circuit("eraser")

    return list(ops)


# ============================================================
# 5. STATE PREPARATION MATRIX
# ============================================================

@lru_cache(maxsize=None)
def householder(prep):
    v = np.zeros(
        2**module.num_weight,
        dtype=complex
    )

    v[:len(prep)] = prep

    e0 = np.zeros_like(v)
    e0[0] = 1

    delta = e0 - v

    if np.linalg.norm(delta) < 1e-14:
        return np.eye(
            len(v),
            dtype=complex
        )

    u = delta / np.linalg.norm(delta)

    return (
        np.eye(len(v), dtype=complex)
        - 2 * np.outer(u, u.conj())
    )


# ============================================================
# 6. COMPUTATIONAL BASIS INDEXING
# ============================================================

@lru_cache(maxsize=None)
def index_pairs(N, target, controls, vals):
    ids = np.arange(N, dtype=np.int32)

    qcount = int(np.log2(N))

    mask = 1 << (qcount - 1 - target)

    selected = (ids & mask) == 0

    for c, v in zip(controls, vals):
        cmask = 1 << (qcount - 1 - c)

        selected &= (
            ((ids & cmask) != 0) == bool(v)
        )

    idx0 = ids[selected]

    return idx0, idx0 | mask


# ============================================================
# 7. APPLY LOGICAL QUANTUM GATES
# ============================================================

def apply(state, op, N):

    if op.typ == "PREP":
        assert op.wires == tuple(
            range(module.num_weight)
        )

        U = householder(op.prep)

        state[:] = (
            U @ state.reshape(len(U), -1)
        ).ravel()

        return

    assert len(op.wires) == 1, op

    target = op.wires[0]

    i, j = index_pairs(
        N,
        target,
        op.controls,
        op.values
    )

    if op.typ == "X":
        b = state[i].copy()
        state[i] = state[j]
        state[j] = b

    elif op.typ == "Z":
        state[j] *= -1

    elif op.typ in ("RY", "H"):

        if op.typ == "H":
            c = 1 / np.sqrt(2)
            s = c
            d = -c

        else:
            c = np.cos(op.theta / 2)
            s = np.sin(op.theta / 2)
            d = c

        a = state[i].copy()
        b = state[j].copy()

        if op.typ == "RY":
            state[i] = c * a - s * b
        else:
            state[i] = c * a + s * b

        state[j] = s * a + d * b

    else:
        raise ValueError(op)


# ============================================================
# 8. PAULI NOISE OPERATIONS
# ============================================================

def apply_pauli(state, wire, ptype, N):

    i, j = index_pairs(
        N,
        wire,
        (),
        ()
    )

    if ptype == "X":
        b = state[i].copy()
        state[i] = state[j]
        state[j] = b

    elif ptype == "Z":
        state[j] *= -1

    elif ptype == "Y":
        a = state[i].copy()
        b = state[j].copy()

        state[i] = -1j * b
        state[j] = 1j * a


# ============================================================
# 9. STATEVECTOR SIMULATION WITH PAULI NOISE
# ============================================================

def simulate(steps, N, prob=0.0, rng=None):

    s = np.zeros(N, dtype=complex)
    s[0] = 1

    if rng is None:
        rng = np.random.default_rng(20261111)

    for op in steps:

        # Apply ideal quantum operation
        apply(s, op, N)

        if prob:

            # Apply stochastic Pauli noise independently
            # to each wire involved in the logical operation.

            wlist = op.wires + op.controls

            flips = rng.random(len(wlist)) < prob

            for wire in np.array(wlist)[flips]:

                ptype = rng.choice(
                    ["X", "Y", "Z"]
                )

                apply_pauli(
                    s,
                    int(wire),
                    ptype,
                    N
                )

    return s


# ============================================================
# 10. EXTRACT OUTPUT DISTRIBUTIONS
# ============================================================

def distribution(state, N, method):

    nW = module.num_weight
    m = module.m
    n = module.n

    if method == "bwsd-min":

        a = (np.abs(state)**2).reshape(
            2**nW,
            2**m,
            2**n,
            2
        )

        result = a.sum(axis=(0, 3))

        return result.sum(axis=0), None

    else:

        p = (np.abs(state)**2).reshape(
            2**nW,
            2**m,
            2**n,
            2
        )

        joint = p.sum(axis=(0, 3))

        output = np.zeros(2**n)

        acc = 0.0

        for zi in range(2**m):

            z = np.array(
                [int(k) for k in f"{zi:0{m}b}"],
                int
            )

            # Accept even-parity outcomes
            if np.sum(z) % 2:
                continue

            # Feedforward correction
            t = module.feedforward_correction(z)

            shift = int(
                "".join(map(str, t)),
                2
            )

            for xi in range(2**n):
                output[xi ^ shift] += joint[zi, xi]

            acc += joint[zi].sum()

        if acc > 1e-13:
            output /= acc

        return output, acc


# ============================================================
# 11. CIRCUIT RESOURCE ANALYSIS
# ============================================================

def describe(name, steps):

    active = set()
    depths = {}

    sizes = Counter()
    types = Counter()

    for op in steps:

        ws = set(
            op.controls + op.wires
        )

        types[op.typ] += 1
        sizes[len(ws)] += 1

        active |= ws

        d = max(
            [depths.get(w, 0) for w in ws],
            default=0
        ) + 1

        for w in ws:
            depths[w] = d

    print(
        f"{name:18s} "
        f"logical gates={len(steps):5d}, "
        f"logical depth={max(depths.values()):5d}, "
        f"wires={len(active):2d}, "
        f"1q={sizes[1]:5d}, "
        f"2q={sizes[2]:5d}, "
        f">2q={sum(z for w, z in sizes.items() if w > 2):5d}, "
        f"max_gate_arity={max(sizes)}"
    )

    print("  types:", dict(types))


# ============================================================
# 12. BUILD CIRCUITS
# ============================================================

circuits = {
    "BWSD-min": (
        "bwsd-min",
        make_ops("bwsd-min")
    ),

    "Grover(0)": (
        "eraser",
        make_ops("eraser", 0)
    ),

    "Grover(1)": (
        "eraser",
        make_ops("eraser", 1)
    ),

    "Grover(4)": (
        "eraser",
        make_ops("eraser", 4)
    )
}


print(
    "SOURCE:",
    args.source,
    "ell=",
    args.ell
)

for name, (_, steps) in circuits.items():
    describe(name, steps)


# ============================================================
# 13. EXACT STATEVECTOR VERIFICATION
# ============================================================

print("\nEXACT STATEVECTOR TESTS")

S = np.array([
    module.score(x)
    for x in __import__("itertools").product(
        (0, 1),
        repeat=module.n
    )
])

best = max(S)

ideal = {}

for name, (method, steps) in circuits.items():

    tic = time.monotonic()

    s = simulate(
        steps,
        2**module.num_qubits
    )

    full = np.abs(
        s.reshape(
            2**module.num_weight,
            2**module.m,
            2**module.n,
            2
        )
    )**2

    unwanted = float(
        full[1:, :, :, :].sum()
        + full[0, :, :, 1].sum()
    )

    assert unwanted < 1e-8, (
        name,
        unwanted
    )

    probs, acc = distribution(
        s,
        len(s),
        method
    )

    opt = float(
        probs[S == best].sum()
    )

    yield_p = opt * (
        acc if acc is not None else 1
    )

    ideal[name] = (
        opt,
        acc,
        yield_p
    )

    print(
        f"{name:18s}: "
        f"P(opt|accepted)={opt:.7f}, "
        f"acceptance={acc if acc is not None else 1:.7f}, "
        f"per-execution optimal yield={yield_p:.7f} "
        f"(t={time.monotonic()-tic:.1f}s)"
    )


# ============================================================
# 14. VERIFY GROVER AMPLIFICATION FORMULA
# ============================================================

p0 = ideal["Grover(0)"][1]

for r in [1, 4]:

    target = np.sin(
        (2*r + 1)
        * np.arcsin(np.sqrt(p0))
    )**2

    actual = ideal[f"Grover({r})"][1]

    print(
        f"Grover {r} "
        f"predicted {target:.10f}, "
        f"simulated {actual:.10f}, "
        f"absolute difference {abs(actual-target):.3e}"
    )

    assert abs(actual - target) < 1e-8

    assert abs(
        ideal[f"Grover({r})"][0]
        - ideal["Grover(0)"][0]
    ) < 1e-8


# ============================================================
# 15. STOCHASTIC PAULI NOISE EXPERIMENT
# ============================================================

rng = np.random.default_rng(args.seed)

rows = []

print(
    "\nSTOCHASTIC LOGICAL-OP PAULI NOISE "
    "(illustrative, NOT physical-gate noise)"
)

for eps in [
    float(x)
    for x in args.error_rates.split(",")
]:

    for name in [
        "BWSD-min",
        "Grover(1)",
        "Grover(4)"
    ]:

        method, steps = circuits[name]

        results = []

        num_runs = (
            1 if eps == 0
            else args.trajectories
        )

        for _ in range(num_runs):

            s = simulate(
                steps,
                2**module.num_qubits,
                eps,
                rng
            )

            probs, acc = distribution(
                s,
                len(s),
                method
            )

            opt = float(
                probs[S == best].sum()
            )

            acceptance = (
                acc if acc is not None else 1
            )

            results.append((
                acceptance,
                opt * acceptance
            ))

        values = np.asarray(results)

        mean_acceptance = values[:, 0].mean()
        mean_yield = values[:, 1].mean()

        yield_sd = (
            values[:, 1].std(ddof=1)
            if len(values) > 1
            else 0
        )

        print(
            f"eps={eps:8.5f} "
            f"{name:18s} "
            f"acc={mean_acceptance:.4%} "
            f"yield={mean_yield:.4%} "
            f"yield_sd={yield_sd:.4%} "
            f"n={len(values)}"
        )

        rows.append({
            "error_rate": eps,
            "method": name,
            "trajectories": len(values),
            "acceptance": mean_acceptance,
            "optimal_yield": mean_yield,
            "yield_std": yield_sd
        })

methods = [
    "BWSD-min",
    "Grover(1)",
    "Grover(4)"
]

colors = {
    "BWSD-min": "royalblue",
    "Grover(1)": "orange",
    "Grover(4)": "green"
}

fig, axes = plt.subplots(
    2, 1,
    figsize=(11, 9),
    sharex=True
)

ax = axes[0]

for method in methods:

    data = sorted(
        [r for r in rows if r["method"] == method],
        key=lambda r: r["error_rate"]
    )

    noise = np.array([
        r["error_rate"] * 100 for r in data
    ])

    optimal_yield = np.array([
        r["optimal_yield"] * 100 for r in data
    ])

    # Approximate 95% Monte Carlo confidence interval
    errors = np.array([
        1.96 * r["yield_std"]
        / np.sqrt(r["trajectories"]) * 100
        for r in data
    ])

    ax.errorbar(
        noise,
        optimal_yield,
        yerr=errors,
        marker="o",
        capsize=4,
        linewidth=2,
        markersize=6,
        label=method,
        color=colors[method]
    )

ax.set_title(
    "DQI Noise Robustness — Optimal Solution Yield",
    fontsize=14
)

ax.set_ylabel("Optimal Solution Yield (%)")
ax.set_ylim(0, 100)
ax.grid(True, linestyle="--", alpha=0.3)
ax.legend()

ax = axes[1]

for method in methods:

    data = sorted(
        [r for r in rows if r["method"] == method],
        key=lambda r: r["error_rate"]
    )

    noise = np.array([
        r["error_rate"] * 100 for r in data
    ])

    acceptance = np.array([
        r["acceptance"] * 100 for r in data
    ])

    ax.plot(
        noise,
        acceptance,
        marker="o",
        linewidth=2,
        markersize=6,
        label=method,
        color=colors[method]
    )

ax.set_title(
    "DQI Noise Robustness — Acceptance Rate",
    fontsize=14
)

ax.set_xlabel("Pauli Error Probability per Wire per Logical Operation (%)")
ax.set_ylabel("Acceptance Rate (%)")
ax.set_ylim(0, 105)

ax.grid(True, linestyle="--", alpha=0.3)
ax.legend()

fig.suptitle(
    f"DQI Max-XORSAT | ell={args.ell} | "
    f"Noise Robustness Comparison",
    fontsize=16,
    fontweight="bold"
)

plt.tight_layout()
plt.savefig(
    f"results/max-xorsat/noise/dqi_noise_ell{args.ell}.png",
    dpi=300,
    bbox_inches="tight"
)
plt.show()
