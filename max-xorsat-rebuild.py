"""Compare collision workarounds for PennyLane DQI Max-XORSAT.

Install: python -m pip install pennylane numpy matplotlib
Run:     python dqi_collision_workarounds.py --ell 3
Optional:python dqi_collision_workarounds.py --ell 3 --shots 100000

The eraser applies Hadamards to the error register. The baseline accepts
only |0...0>. This extended experiment accepts ALL even-parity outcomes
and corrects the final assignment using classical feedforward. An optional
Grover amplitude amplification circuit boosts the even-parity event.
Neither method is guaranteed to outperform a standard DQI decoder.
"""
from __future__ import annotations

import argparse
from itertools import combinations, product
from math import ceil, comb, log2

import matplotlib.pyplot as plt
import numpy as np
import pennylane as qml

parser = argparse.ArgumentParser()
parser.add_argument('--ell', type=int, default=3, help='Polynomial degree, 0..5')
parser.add_argument('--shots', type=int, default=0, help='Optional finite-shot eraser test; e.g. 100000')
parser.add_argument('--save', type=str, default='', help='Optional path to save plot (e.g. comparison.png)')
parser.add_argument('--grover', type=int, default=0, help='Optional Grover iterations boosting even-parity feedforward success')
args = parser.parse_args()

B = np.array([
    [1, 0, 0, 0],
    [1, 1, 0, 0],
    [0, 1, 1, 0],
    [0, 0, 1, 1],
    [0, 0, 0, 1],
], dtype=int)
v = np.array([1, 0, 1, 0, 1], dtype=int)
m, n = B.shape
ell = args.ell
if not (0 <= ell <= m):
    parser.error(f'--ell must be between 0 and {m}')
if args.shots < 0:
    parser.error('--shots must be >= 0')
if args.grover < 0:
    parser.error('--grover must be >= 0')

num_weight = max(1, ceil(log2(ell + 1)))  # 1 idle qubit at ell=0
weight_wires = list(range(num_weight))
error_wires = list(range(num_weight, num_weight + m))
syndrome_wires = list(range(num_weight + m, num_weight + m + n))
match_wire = num_weight + m + n
num_qubits = match_wire + 1


def score(x):
    x = np.asarray(x, dtype=int)
    return int(np.count_nonzero((B @ x) % 2 == v))


def optimal_weights():
    A = np.zeros((ell + 1, ell + 1))
    for k in range(1, ell + 1):
        z = np.sqrt(k * (m - k + 1))
        A[k - 1, k] = z
        A[k, k - 1] = z
    _, vecs = np.linalg.eigh(A)
    w = vecs[:, -1]
    return w if w[0] >= 0 else -w


weights = optimal_weights()


def split_cycle_shift(active_size, excitations):
    last = error_wires[active_size - 1]
    qml.CNOT(wires=[last - 1, last])
    theta = 2 * np.arccos(np.sqrt(1 / active_size))
    qml.CRY(theta, wires=[last, last - 1])
    qml.CNOT(wires=[last - 1, last])
    for r in range(2, excitations + 1):
        qml.CNOT(wires=[last - r, last])
        theta = 2 * np.arccos(np.sqrt(r / active_size))
        qml.ctrl(qml.RY, control=[last, last - r + 1])(
            theta, wires=last - r
        )
        qml.CNOT(wires=[last - r, last])


def prepare_dicke(k):
    for wire in error_wires[m - k:]:
        qml.PauliX(wires=wire)
    for size in range(m, k, -1):
        split_cycle_shift(size, k)
    for size in range(k, 1, -1):
        split_cycle_shift(size, size - 1)


def prepare_weight_error():
    qml.StatePrep(weights, wires=weight_wires, pad_with=0.0)
    for k in range(1, ell + 1):
        k_bits = [int(bit) for bit in format(k, f'0{num_weight}b')]
        qml.ctrl(prepare_dicke, control=weight_wires, control_values=k_bits)(k)


def erase_weight_register():
    for k in range(1, ell + 1):
        k_bits = [int(bit) for bit in format(k, f'0{num_weight}b')]
        for chosen in combinations(range(m), k):
            error_bits = [int(i in chosen) for i in range(m)]
            for j, bit in enumerate(k_bits):
                if bit:
                    qml.ctrl(
                        qml.PauliX, control=error_wires, control_values=error_bits
                    )(wires=weight_wires[j])


def encode_constraint_phases():
    for i, bit in enumerate(v):
        if bit:
            qml.PauliZ(wires=error_wires[i])


def compute_syndrome():
    for j in range(m):
        for i in range(n):
            if B[j, i]:
                qml.CNOT(wires=[error_wires[j], syndrome_wires[i]])


def build_bwsd_candidates(max_weight):
    groups = {}
    for k in range(max_weight + 1):
        for chosen in combinations(range(m), k):
            error = tuple(int(j in chosen) for j in range(m))
            syndrome = tuple(int(x) for x in ((B.T @ np.asarray(error)) % 2))
            groups.setdefault(syndrome, []).append(error)
    return groups


def select_corrections(policy='min'):
    groups = build_bwsd_candidates(ell)
    key = lambda e: (sum(e), e)
    if policy not in ('min', 'max'):
        raise ValueError(f'Unknown correction policy: {policy}')
    chooser = min if policy == 'min' else max
    return {syndrome: chooser(errors, key=key) for syndrome, errors in groups.items()}


def apply_bwsd(policy):
    for syndrome, selected in sorted(select_corrections(policy).items()):
        if not any(selected):
            continue
        qml.ctrl(
            qml.PauliX, control=syndrome_wires, control_values=syndrome
        )(wires=match_wire)
        for j, bit in enumerate(selected):
            if bit:
                qml.CNOT(wires=[match_wire, error_wires[j]])
        qml.ctrl(
            qml.PauliX, control=syndrome_wires, control_values=syndrome
        )(wires=match_wire)


def dqi_circuit(method):
    prepare_weight_error()
    erase_weight_register()
    encode_constraint_phases()
    compute_syndrome()

    if method == 'eraser':
        # Change to the X basis; selecting error=0 aggregates all error
        # branches with the same syndrome, *with their signed amplitudes*.
        for wire in error_wires:
            qml.Hadamard(wires=wire)
    elif method == 'bwsd-min':
        apply_bwsd('min')
    elif method == 'bwsd-max':
        apply_bwsd('max')
    else:
        raise ValueError(f'Unknown method: {method}')

    for wire in syndrome_wires:
        qml.Hadamard(wires=wire)


exact_device = qml.device('default.qubit', wires=num_qubits)


@qml.qnode(exact_device)
def joint_probabilities(method):
    dqi_circuit(method)
    return qml.probs(wires=error_wires + syndrome_wires)


def postselect_error_zero(method):
    """Return conditional P(x) and *physical* acceptance probability."""
    joint = np.asarray(joint_probabilities(method), dtype=float).reshape(2**m, 2**n)
    selected = joint[0]  # error register = |00000> after BWSD or H
    acceptance = float(np.sum(selected))
    if acceptance < 1e-12:
        return None, acceptance
    return selected / acceptance, acceptance


def unconditional_probabilities(method):
    joint = np.asarray(joint_probabilities(method), dtype=float).reshape(2**m, 2**n)
    return joint.sum(axis=0)


def classical_caa_reference():
    """Exponential reference, NOT an efficient quantum algorithm.

    Accumulate signed amplitudes at each common syndrome, Hadamard transform,
    and normalize the resulting polynomial distribution. Also return the
    expected quantum-eraser acceptance probability (i.e. /2^m).
    """
    a_s = np.zeros(2**n, dtype=float)
    for y in product((0, 1), repeat=m):
        k = sum(y)
        if k > ell:
            continue
        error = np.asarray(y, dtype=int)
        syndrome = (B.T @ error) % 2
        idx = int(''.join(map(str, syndrome)), 2)
        signed = (-1) ** int(v @ error)
        a_s[idx] += signed * weights[k] / np.sqrt(comb(m, k))

    # Hadamard transform of the aggregate: sum_s (-1)^{x.s} a_s / sqrt(2^n)
    a_x = np.zeros(2**n, dtype=float)
    for x in product((0, 1), repeat=n):
        idx = int(''.join(map(str, x)), 2)
        a_x[idx] = sum(
            a_s[s_idx] * (-1) ** (int(''.join(map(str, x)), 2) & s_idx).bit_count()
            for s_idx in range(2**n)
        ) / np.sqrt(2**n)
    raw = a_x**2
    norm = float(np.sum(raw))
    success = norm / (2**m)
    if norm < 1e-12:
        return None, success
    return raw / norm, success



def reflect_initial_zero():
    """I - 2|0...0><0...0| on all *active* qubits.

    The final match ancilla is unused by the eraser and may be omitted.
    """
    active = weight_wires + error_wires + syndrome_wires
    for wire in active:
        qml.PauliX(wires=wire)
    qml.ctrl(qml.PauliZ, control=active[:-1])(wires=active[-1])
    for wire in active:
        qml.PauliX(wires=wire)


@qml.qnode(exact_device)
def amplified_eraser_joint(grover_steps):
    """Grover amplitude amplification of even-parity error outcomes.

    The even-parity phase oracle is (-Z)^tensor(m), equal up to global
    phase to applying Z to each error qubit. Requires full A and A^dagger.
    """
    dqi_circuit('eraser')
    for _ in range(grover_steps):
        # Oracle: phase distinguishes even and odd parity of error register
        for wire in error_wires:
            qml.PauliZ(wires=wire)
        qml.adjoint(dqi_circuit)('eraser')
        reflect_initial_zero()
        dqi_circuit('eraser')
    return qml.probs(wires=error_wires + syndrome_wires)


def feedforward_correction(error_bits):
    """For this exact 5x4 B: solve B @ t = z over GF(2).

    On even-parity z, B has a unique preimage t. The prefix-XOR
    identity below is specific to this B, not a generic GF(2) solver.
    """
    z = np.asarray(error_bits, dtype=int)
    if z.shape != (m,) or int(z.sum()) % 2:
        raise ValueError("Expected one even-parity length-m error-register outcome")
    t = np.bitwise_xor.accumulate(z[:n])
    np.testing.assert_array_equal((B @ t) % 2, z)
    return t


def feedforward_postselect_exact(joint):
    """Accept *all* z in Im(B), correct measured x to x XOR t(z)."""
    joint = np.asarray(joint, dtype=float).reshape(2**m, 2**n)
    raw = np.zeros(2**n, dtype=float)
    for zi in range(2**m):
        z = np.array([int(q) for q in f'{zi:0{m}b}'], dtype=int)
        if int(z.sum()) % 2:
            continue
        t = feedforward_correction(z)
        shift = int(''.join(map(str, t)), 2)
        for xi in range(2**n):
            raw[xi ^ shift] += joint[zi, xi]
    acceptance = float(raw.sum())
    if acceptance < 1e-12:
        return None, acceptance
    return raw / acceptance, acceptance


def feedforward_postselect_samples(samples):
    """Convert joint Z-basis samples taken after the eraser Hadamards.

    Sample row ordering: error-register 5 bits, then syndrome 4 bits.
    """
    samples = np.asarray(samples, dtype=int).reshape(-1, m + n)
    zs, xs = samples[:, :m], samples[:, m:]
    mask = (zs.sum(axis=1) % 2 == 0)
    zs, xs = zs[mask], xs[mask]
    if len(xs) == 0:
        return None, 0
    ts = np.bitwise_xor.accumulate(zs[:, :n], axis=1)
    corrected_x = np.bitwise_xor(xs, ts)
    indices = corrected_x @ (2 ** np.arange(n - 1, -1, -1))
    probs = np.bincount(indices, minlength=2**n) / len(indices)
    return probs, int(len(indices))

def analyze_method(label, probs, acceptance, scores, best):
    if probs is None:
        print(f'{label:26s} | acceptance={acceptance:.8%} | conditional distribution undefined')
        return
    assert np.isclose(np.sum(probs), 1.0, atol=1e-9)
    p_opt = float(sum(p for p, s in zip(probs, scores) if s == best))
    avg = float(np.dot(probs, scores))
    print(f'{label:26s} | P(opt)={p_opt:8.4%} | '
          f'mean score={avg:.5f} | acceptance={acceptance:9.5%}')


def main():
    bitstrings = [''.join(map(str, x)) for x in product((0, 1), repeat=n)]
    scores = [score(x) for x in product((0, 1), repeat=n)]
    best = max(scores)
    uniform = np.ones(2**n) / 2**n

    groups = build_bwsd_candidates(ell)
    collisions = {s: ys for s, ys in groups.items() if len(ys) > 1}
    print(f'DQI Max-XORSAT | ell={ell}, m={m}, n={n}, qubits={num_qubits}')
    print(f'Prepared errors: {sum(comb(m, k) for k in range(ell + 1))}')
    print(f'Distinct syndromes: {len(groups)}, colliding groups: {len(collisions)}')
    print('Optimal weights:', np.round(weights, 8))
    print(f'Best score = {best}/{m}\n')

    p_min = unconditional_probabilities('bwsd-min')
    p_max = unconditional_probabilities('bwsd-max')
    p_min_post, acc_min = postselect_error_zero('bwsd-min')
    p_eraser, acc_eraser = postselect_error_zero('eraser')
    p_caa, acc_caa = classical_caa_reference()
    joint_eraser = np.asarray(joint_probabilities('eraser'), dtype=float)
    p_ff, acc_ff = feedforward_postselect_exact(joint_eraser)

    np.testing.assert_allclose(p_min.sum(), 1.0, atol=1e-9)
    np.testing.assert_allclose(p_max.sum(), 1.0, atol=1e-9)
    np.testing.assert_allclose(acc_eraser, acc_caa, atol=1e-10)
    if p_caa is not None:
        np.testing.assert_allclose(p_eraser, p_caa, atol=1e-9)
        np.testing.assert_allclose(p_ff, p_eraser, atol=1e-9)
    np.testing.assert_allclose(acc_ff, (2**n) * acc_eraser, atol=1e-9)
    print('Quantum eraser vs classical CAA: PASS\n')

    analyze_method('Uniform', uniform, 1.0, scores, best)
    analyze_method('BWSD-min (all shots)', p_min, 1.0, scores, best)
    analyze_method('BWSD-max (all shots)', p_max, 1.0, scores, best)
    analyze_method('BWSD-min (postselected)', p_min_post, acc_min, scores, best)
    analyze_method('Quantum eraser (postsel.)', p_eraser, acc_eraser, scores, best)
    analyze_method('Eraser + feedforward', p_ff, acc_ff, scores, best)

    acc_amp = None
    if args.grover:
        joint_amp = amplified_eraser_joint(args.grover)
        p_amp, acc_amp = feedforward_postselect_exact(joint_amp)
        predicted = np.sin((2 * args.grover + 1) * np.arcsin(np.sqrt(acc_ff))) ** 2
        np.testing.assert_allclose(acc_amp, predicted, atol=1e-8)
        np.testing.assert_allclose(p_amp, p_ff, atol=1e-8)
        analyze_method(f'Eraser + FF + Grover({args.grover})', p_amp, acc_amp, scores, best)
    analyze_method('Classical CAA reference', p_caa, acc_caa, scores, best)

    if args.shots > 0:
        shot_device = qml.device('default.qubit', wires=num_qubits, shots=args.shots)

        @qml.qnode(shot_device)
        def sample_eraser_joint():
            dqi_circuit('eraser')
            return qml.sample(wires=error_wires + syndrome_wires)

        samples = np.asarray(sample_eraser_joint(), dtype=int).reshape(-1, m + n)
        keep = samples[np.all(samples[:, :m] == 0, axis=1), m:]
        print(f'\nFinite shots quantum eraser: accepted {len(keep)} / {args.shots} '
              f'= {len(keep) / args.shots:.4%}')
        if len(keep):
            sampled_x = keep @ (2 ** np.arange(n - 1, -1, -1))
            sampled_probs = np.bincount(sampled_x, minlength=2**n) / len(keep)
            analyze_method('Quantum eraser (sampled)', sampled_probs,
                           len(keep) / args.shots, scores, best)
        sampled_ff_probs, n_ff = feedforward_postselect_samples(samples)
        print(f'Feedforward: accepted {n_ff} / {args.shots} = {n_ff / args.shots:.4%}')
        analyze_method('Eraser + feedforward (sampled)', sampled_ff_probs,
                       n_ff / args.shots, scores, best)

    # Plot conditional distributions explicitly; compare acceptance separately.
    curves = [
        ('Uniform', uniform, 'gray'),
        ('BWSD-min', p_min, '#4472C4'),
        ('BWSD-max', p_max, '#ED7D31'),
        ('BWSD-min + postsel.', p_min_post, '#70AD47'),
        ('Quantum eraser + postsel.', p_eraser, '#7030A0'),
    ]
    curves = [(a, b, c) for a, b, c in curves if b is not None]
    positions = np.arange(2**n)
    width = 0.78 / len(curves)
    fig, (ax, ax_accept) = plt.subplots(2, 1, figsize=(14, 8),
                                       gridspec_kw={'height_ratios': [3, 1]})
    for i, (name, probs, color) in enumerate(curves):
        ax.bar(positions - 0.39 + (i + 0.5)*width, probs, width=width,
               label=name, color=color, alpha=0.85)
    ax.set_xticks(positions)
    ax.set_xticklabels(bitstrings, rotation=45)
    for tick, s in zip(ax.get_xticklabels(), scores):
        if s == best:
            tick.set_color('green')
            tick.set_fontweight('bold')
    ax.set_ylabel('P(assignment | accepted, when applicable)')
    ax.set_title(f'DQI collision workarounds — ell={ell} (feedforward has same purple distribution)')
    ax.grid(axis='y', alpha=0.25, linestyle='--')
    ax.legend(ncol=2, fontsize=9)

    accept_names = ['BWSD-min postsel.', 'Original eraser', 'Eraser + feedforward']
    accept_probs = [acc_min, acc_eraser, acc_ff]
    accept_colors = ['#70AD47', '#7030A0', '#C39BD3']
    if acc_amp is not None:
        accept_names.append(f'FF + Grover({args.grover})')
        accept_probs.append(acc_amp)
        accept_colors.append('#2266AA')
    ax_accept.bar(accept_names, accept_probs, color=accept_colors)
    ax_accept.set_ylabel('Acceptance rate')
    ax_accept.set_ylim(0, 1.1)
    for i, a in enumerate(accept_probs):
        ax_accept.text(i, a + 0.02, f'{a:.2%}', ha='center', fontsize=10)
    ax_accept.grid(axis='y', alpha=0.25)
    fig.tight_layout()
    if args.save:
        fig.savefig(args.save, dpi=180)
        print(f'\nSaved plot: {args.save}')
    plt.show()


if __name__ == '__main__':
    main()
