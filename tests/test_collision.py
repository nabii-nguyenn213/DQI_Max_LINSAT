from ast import arg
from pathlib import Path
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

from src.utils import check_collisions

ROOT = Path(__file__).resolve().parent.parent

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

    args = parser.parse_args()

    if args.shots <= 0:
        parser.error("--shots must be positive")

    if args.degree < 0:
        parser.error("--degree must be nonnegative")

    q, B, allowed, weights = load_instance(
        args.dataset
    )

    check_collisions(B, q, degree=args.degree)

if __name__ == "__main__": 
    main()
