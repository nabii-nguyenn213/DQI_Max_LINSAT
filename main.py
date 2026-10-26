from pathlib import Path
from itertools import product
from collections import Counter
from time import perf_counter
import argparse

import numpy as np
import pandas as pd

from src.formulation import (
    get_matrices,
    find_finite_field,
    eligible_lecturers,
    generate_conflicts,
    formulate_max_linsat,
)

from src.dqi import DQI
from src.visualize import DQIVisualizer


ROOT = Path(__file__).resolve().parent


# ==========================================
# 1. Load dataset and formulate Max-LINSAT
# ==========================================

def load_instance(dataset="tiny"):

    data_dir = ROOT / "datasets" / dataset

    if not data_dir.is_dir():
        raise FileNotFoundError(
            f"Dataset directory not found: {data_dir}"
        )

    skills, slots, _, tasks, subject_ids, slot_ids = (
        get_matrices(
            skillPath=str(data_dir / "InstructorSkill.csv"),
            slotPath=str(data_dir / "InstructorSlot.csv"),
            conflictPath=str(data_dir / "SlotConflict.csv"),
            taskPath=str(data_dir / "Task.csv"),
        )
    )

    # Read slot-conflict adjacency matrix
    slot_ids = [str(s) for s in slot_ids]

    conflict_df = pd.read_csv(
        data_dir / "SlotConflict.csv",
        index_col=0
    )

    conflict_df.index = conflict_df.index.map(str)
    conflict_df.columns = conflict_df.columns.map(str)

    slot_matrix = conflict_df.loc[
        slot_ids, slot_ids
    ].to_numpy(dtype=int)

    # Determine finite field
    q = find_finite_field(skills.shape[0])

    # Generate eligibility constraints
    eligible = eligible_lecturers(
        skills, slots, tasks, subject_ids, slot_ids
    )

    # Generate conflict constraints
    conflicts = generate_conflicts(
        tasks, slot_matrix, slot_ids
    )

    # Construct Max-LINSAT
    B, allowed, weights, labels = formulate_max_linsat(
        q, eligible, conflicts
    )

    print("\n=== Max-LINSAT formulation ===")
    print("Dataset:", dataset)
    print("Field:", q)
    print("Variables:", B.shape[1])
    print("Constraints:", B.shape[0])
    print("Eligible:", eligible)
    print("Conflicts:", conflicts)
    print("\nB:\n", B)
    print("\nAllowed:\n", allowed.astype(int))

    # Exact validation only for original tiny dataset
    if dataset == "tiny":

        assert conflicts == [
            (0, 1), (0, 2), (1, 2)
        ]

        assert B.shape == (6, 3)

    return q, B, allowed, weights


# ==========================================
# 2. Run DQI with selected decoder
# ==========================================

def run_decoder(dqi, name, shots=2000, seed=42):

    print(f"\n=== Running {name.upper()} ===")

    # Execute PennyLane DQI.
    # The DQI class creates the selected decoder,
    # optimizes its polynomial coefficients,
    # and performs postselected decoding.

    start = perf_counter()

    probabilities = dqi.output_probabilities(
        decoder_type=name
    )

    runtime = perf_counter() - start

    # IMPORTANT:
    # Save this immediately before another
    # decoder changes the DQI class state.
    success_probability = (
        dqi.last_success_probability
    )

    if success_probability is None:
        raise RuntimeError(
            "DQI did not report decoding success"
        )

    # --------------------------------------
    # 3. Decoder-specific classical reference
    # --------------------------------------

    decoder = dqi.create_decoder(name)

    # Uses the coefficients optimized for
    # this decoder by output_probabilities().
    reference = dqi.reference_probabilities(
        decoder
    )

    # --------------------------------------
    # 4. Enumerate assignments
    # --------------------------------------

    assignments = list(
        product(range(dqi.q), repeat=dqi.n)
    )

    p = np.array([
        probabilities.get(y, 0.0)
        for y in assignments
    ], dtype=float)

    if np.any(p < -1e-10):
        raise ValueError(
            f"{name}: negative probabilities"
        )

    if not np.isclose(p.sum(), 1.0, atol=1e-8):
        raise ValueError(
            f"{name}: probabilities not normalized"
        )

    p = np.clip(p, 0.0, None)
    p /= p.sum()

    probabilities = dict(
        zip(assignments, p.tolist())
    )

    # --------------------------------------
    # 5. Verify PennyLane against reference
    # --------------------------------------

    difference = max(
        abs(probabilities[y] - reference[y])
        for y in assignments
    )

    if difference > 1e-6:
        raise AssertionError(
            f"{name} differs from reference: "
            f"{difference}"
        )

    # --------------------------------------
    # 6. Finite-shot sampling
    # --------------------------------------

    rng = np.random.default_rng(seed)

    samples = rng.choice(
        len(assignments),
        size=shots,
        p=p
    )

    counts = Counter(
        assignments[int(i)]
        for i in samples
    )

    # --------------------------------------
    # 7. Calculate scores
    # --------------------------------------

    scores = np.array([
        dqi.objective(y)
        for y in assignments
    ])

    best_score = float(scores.max())

    exact_score = float(
        np.dot(p, scores)
    )

    sampled_score = float(
        np.mean(scores[samples])
    )

    optimal_probability = float(
        p[np.isclose(scores, best_score)].sum()
    )

    # --------------------------------------
    # 8. Store results
    # --------------------------------------

    result = {
        "probabilities": probabilities,
        "reference": reference,
        "counts": dict(counts),
        "exact_score": exact_score,
        "sampled_score": sampled_score,
        "optimal_probability": optimal_probability,
        "success_probability": float(
            success_probability
        ),
        "runtime_seconds": runtime,
        "reference_difference": difference,
    }

    # --------------------------------------
    # 9. Print results
    # --------------------------------------

    print(f"Expected score: {exact_score:.6f}")
    print(f"Sampled score: {sampled_score:.6f}")

    print(
        f"Optimal probability: "
        f"{optimal_probability:.4%}"
    )

    print(
        f"Postselection success: "
        f"{success_probability:.4%}"
    )

    print(
        f"DQI execution time: {runtime:.4f} s"
    )

    print(
        f"Reference difference: {difference:.3e}"
    )

    print("\nTop assignments:")

    for assignment, count in counts.most_common(5):
        print(f"  {assignment}: {count}")

    return result


# ==========================================
# 3. Compare both decoders
# ==========================================

def compare_decoders(results):

    lut = results["lut"]
    bwsd = results["bwsd"]

    assignments = lut["probabilities"].keys()

    difference = max(
        abs(
            lut["probabilities"][y]
            - bwsd["probabilities"][y]
        )
        for y in assignments
    )

    print("\n=== LUT vs BWSD ===")

    print(
        "Maximum probability difference:",
        f"{difference:.8f}"
    )

    print(
        "Expected score difference:",
        f"{abs(lut['exact_score'] - bwsd['exact_score']):.6f}"
    )

    print(
        "LUT success:",
        f"{lut['success_probability']:.2%}"
    )

    print(
        "BWSD success:",
        f"{bwsd['success_probability']:.2%}"
    )

    print(
        "LUT execution:",
        f"{lut['runtime_seconds']:.4f} s"
    )

    print(
        "BWSD execution:",
        f"{bwsd['runtime_seconds']:.4f} s"
    )


# ==========================================
# 4. Main experiment
# ==========================================

def main():

    parser = argparse.ArgumentParser(
        description="DQI LUT vs BWSD comparison"
    )

    parser.add_argument(
        "--dataset",
        type=str,
        default="tiny"
    )

    parser.add_argument(
        "--degree",
        type=int,
        default=1
    )

    parser.add_argument(
        "--shots",
        type=int,
        default=2000
    )

    parser.add_argument(
        "--no-show",
        action="store_true"
    )

    args = parser.parse_args()

    if args.shots <= 0:
        parser.error("--shots must be positive")

    if args.degree < 0:
        parser.error("--degree must be nonnegative")

    # Step 1: Formulate
    q, B, allowed, weights = load_instance(
        args.dataset
    )

    # Step 2: Initialize DQI
    dqi = DQI(
        B=B,
        allowed=allowed,
        weights=weights,
        q=q,
        degree=args.degree
    )

    # Step 3: Run both decoders
    results = {}

    for name, seed in [
        ("lut", 42),
        ("bwsd", 43)
    ]:

        results[name] = run_decoder(
            dqi=dqi,
            name=name,
            shots=args.shots,
            seed=seed
        )

    # Step 4: Compare results
    compare_decoders(results)

    # Step 5: Visualize
    visualizer = DQIVisualizer(dqi)

    output_path = visualizer.plot_comparison(
        results=results,
        shots=args.shots,
        output_dir=ROOT / "results",
        show=not args.no_show
    )

    print("\nFigure saved to:", output_path)


if __name__ == "__main__":
    main()
