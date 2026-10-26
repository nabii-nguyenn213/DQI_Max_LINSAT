from pathlib import Path
from itertools import product

import numpy as np
import matplotlib.pyplot as plt


class DQIVisualizer:

    def __init__(self, dqi):
        self.dqi = dqi

    def plot_comparison(
        self,
        results,
        reference=None,
        shots=2000,
        output_dir="results",
        show=True
    ):
        """
        Visualize LUT vs BWSD DQI results.

        Each decoder result should contain:
            probabilities: {assignment: probability}
            counts: {assignment: sample_count}
            reference: {assignment: reference_probability}
            success_probability: float

        The reference parameter is retained for
        backward compatibility with degree 1.
        """

        output_dir = Path(output_dir)
        output_dir.mkdir(
            parents=True,
            exist_ok=True
        )

        decoders = ("lut", "bwsd")

        # -------------------------------------
        # 1. Enumerate all assignments
        # -------------------------------------

        assignments = list(
            product(
                range(self.dqi.q),
                repeat=self.dqi.n
            )
        )

        # Sort from highest score to lowest
        assignments.sort(
            key=lambda y: (
                -self.dqi.objective(y), y
            )
        )

        x = np.arange(len(assignments))

        labels = [
            "".join(map(str, y))
            if self.dqi.q <= 10 else str(y)
            for y in assignments
        ]

        scores = np.array([
            self.dqi.objective(y)
            for y in assignments
        ])

        # -------------------------------------
        # 2. Load decoder-specific results
        # -------------------------------------

        probabilities = {}
        references = {}
        counts = {}
        success = {}
        expected_scores = {}

        for name in decoders:

            if name not in results:
                raise ValueError(
                    f"Missing results for {name}"
                )

            result = results[name]

            probabilities[name] = np.array([
                result["probabilities"].get(y, 0.0)
                for y in assignments
            ], dtype=float)

            if not np.isclose(
                probabilities[name].sum(),
                1.0,
                atol=1e-6
            ):
                raise ValueError(
                    f"{name} probabilities are not normalized"
                )

            # Use each decoder's own reference
            decoder_reference = result.get("reference")

            if decoder_reference is None:

                if self.dqi.degree > 1:
                    raise ValueError(
                        f"Missing postselected reference "
                        f"for {name}. Store it in "
                        f"results['{name}']['reference']."
                    )

                decoder_reference = reference

            if decoder_reference is not None:

                references[name] = np.array([
                    decoder_reference.get(y, 0.0)
                    for y in assignments
                ], dtype=float)

            else:
                references[name] = None

            counts[name] = np.array([
                result["counts"].get(y, 0)
                for y in assignments
            ], dtype=int)

            success[name] = result.get(
                "success_probability"
            )

            expected_scores[name] = float(
                np.dot(probabilities[name], scores)
            )

        # Require success rates for postselected DQI
        if self.dqi.degree > 1:
            for name in decoders:
                if success[name] is None:
                    raise ValueError(
                        f"Missing success_probability "
                        f"for {name}"
                    )

        has_success = all(
            success[name] is not None
            for name in decoders
        )

        # -------------------------------------
        # 3. Create figure
        # -------------------------------------

        nplots = 4 if has_success else 3

        height_ratios = (
            [1.6, 1.0, 1.2, 0.7]
            if has_success
            else [1.6, 1.0, 1.2]
        )

        fig, axes = plt.subplots(
            nplots,
            1,
            figsize=(17, 17 if has_success else 13),
            gridspec_kw={
                "height_ratios": height_ratios
            },
            layout="constrained"
        )

        fig.suptitle(
            f"DQI Tiny Dataset | "
            f"q={self.dqi.q} | "
            f"degree={self.dqi.degree} | "
            f"shots={shots}",
            fontsize=16,
            fontweight="bold"
        )

        width = 0.35

        # =====================================
        # Plot 1: Assignment probabilities
        # =====================================

        ax = axes[0]

        ax.bar(
            x - width / 2,
            probabilities["lut"],
            width,
            label="PennyLane + LUT",
            alpha=0.8
        )

        ax.bar(
            x + width / 2,
            probabilities["bwsd"],
            width,
            label="PennyLane + BWSD",
            alpha=0.8
        )

        # Separate reference for each decoder
        markers = {
            "lut": "o",
            "bwsd": "x"
        }

        for name in decoders:

            if references[name] is None:
                continue

            ax.plot(
                x,
                references[name],
                linestyle="None",
                marker=markers[name],
                markersize=5,
                label=f"{name.upper()} reference"
            )

        title = "Assignment probability distribution"

        if self.dqi.degree > 1:
            title += " (conditioned on decoding success)"

        ax.set_title(title)
        ax.set_ylabel("Probability")

        ax.set_xticks(x)
        ax.set_xticklabels(
            labels,
            rotation=70,
            fontsize=8
        )

        ax.legend()
        ax.grid(axis="y", alpha=0.25)

        # =====================================
        # Plot 2: Objective score distribution
        # =====================================

        ax = axes[1]

        unique_scores = np.unique(scores)
        sx = np.arange(len(unique_scores))

        # Uniform baseline
        uniform_mass = [
            np.mean(np.isclose(scores, s))
            for s in unique_scores
        ]

        # Decoder score distributions
        score_mass = {}

        for name in decoders:

            score_mass[name] = [
                probabilities[name][
                    np.isclose(scores, s)
                ].sum()
                for s in unique_scores
            ]

        sw = 0.26

        ax.bar(
            sx - sw,
            uniform_mass,
            sw,
            label="Uniform random"
        )

        ax.bar(
            sx,
            score_mass["lut"],
            sw,
            label=(
                f"LUT (E[S]={expected_scores['lut']:.4f})"
            )
        )

        ax.bar(
            sx + sw,
            score_mass["bwsd"],
            sw,
            label=(
                f"BWSD (E[S]={expected_scores['bwsd']:.4f})"
            )
        )

        ax.set_title(
            "Max-LINSAT objective score distribution"
        )

        ax.set_xlabel("Objective score")
        ax.set_ylabel("Probability mass")

        ax.set_xticks(sx)
        ax.set_xticklabels([
            f"{s:g}" for s in unique_scores
        ])

        ax.legend()
        ax.grid(axis="y", alpha=0.25)

        # =====================================
        # Plot 3: Finite-shot samples
        # =====================================

        ax = axes[2]

        ax.bar(
            x - width / 2,
            counts["lut"],
            width,
            label="LUT samples"
        )

        ax.bar(
            x + width / 2,
            counts["bwsd"],
            width,
            label="BWSD samples"
        )

        if self.dqi.degree > 1:
            ax.set_title(
                "Accepted assignment samples "
                "(postselected)"
            )
        else:
            ax.set_title(
                "Sampled lecturer assignments"
            )

        ax.set_xlabel("Lecturer assignment")
        ax.set_ylabel("Sample count")

        ax.set_xticks(x)
        ax.set_xticklabels(
            labels,
            rotation=70,
            fontsize=8
        )

        ax.legend()
        ax.grid(axis="y", alpha=0.25)

        # =====================================
        # Plot 4: Decoding success probability
        # =====================================

        if has_success:

            ax = axes[3]


            values = np.array([
                float(success["lut"]),
                float(success["bwsd"])
            ], dtype=float)

            # Allow small floating-point numerical errors
            tolerance = 1e-8

            if (
                not np.all(np.isfinite(values))
                or np.any(values < -tolerance)
                or np.any(values > 1.0 + tolerance)
            ):
                raise ValueError(
                    f"Invalid decoding success probability: {values}"
                )

            # Correct only floating-point rounding errors
            values = np.clip(values, 0.0, 1.0)

            bars = ax.bar(
                ["LUT", "BWSD"],
                values,
                width=0.45
            )

            ax.set_title(
                "Decoding / postselection success probability"
            )

            ax.set_ylabel("Success probability")
            ax.set_ylim(0, 1.1)

            for bar, value in zip(bars, values):

                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value + 0.025,
                    f"{value * 100:.2f}%",
                    ha="center",
                    va="bottom",
                    fontweight="bold"
                )

            ax.grid(axis="y", alpha=0.25)

        # -------------------------------------
        # 4. Save figure
        # -------------------------------------

        output = (
            output_dir
            / "tiny_dqi_lut_vs_bwsd.png"
        )

        fig.savefig(
            output,
            dpi=170,
            bbox_inches="tight"
        )

        if show:
            plt.show()

        plt.close(fig)

        return output
