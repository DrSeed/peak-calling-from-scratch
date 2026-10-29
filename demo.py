# Self-contained demo: simulate coverage with embedded peaks, call them, and plot.
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from peak_caller import call_peaks


def simulate_coverage(length=5000, baseline=4.0, seed=42):
    # Build a Poisson background track with a few Gaussian-shaped enrichment peaks.
    rng = np.random.default_rng(seed)
    signal = np.full(length, baseline)
    true_peaks = [(800, 60, 40), (2200, 100, 25), (3600, 40, 60), (4400, 80, 15)]
    positions = np.arange(length)
    for centre, width, height in true_peaks:
        signal = signal + height * np.exp(-0.5 * ((positions - centre) / width) ** 2)
    coverage = rng.poisson(signal).astype(float)
    return coverage, true_peaks


def main():
    os.makedirs("figures", exist_ok=True)
    os.makedirs("results", exist_ok=True)

    coverage, true_peaks = simulate_coverage()
    score_threshold = 6.0
    peaks, scores = call_peaks(
        coverage,
        window=40,
        flank=800,
        score_threshold=score_threshold,
        min_gap=120,
    )

    # Save the summary table.
    if peaks:
        df = pd.DataFrame(peaks)
    else:
        df = pd.DataFrame(columns=["start", "end", "summit", "score", "max_coverage"])
    df.to_csv("results/summary.csv", index=False)

    # Plot coverage, scores, and called peak boundaries.
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 7), sharex=True)

    ax1.plot(coverage, color="#4477aa", lw=0.7, label="Coverage")
    for pk in peaks:
        ax1.axvspan(pk["start"], pk["end"], color="#ee6677", alpha=0.3)
        ax1.axvline(pk["summit"], color="#ee6677", lw=0.8, ls="--")
    ax1.set_ylabel("Read coverage")
    ax1.set_title("Simulated coverage with called peaks (shaded)")
    ax1.legend(loc="upper right")

    ax2.plot(scores, color="#228833", lw=0.7, label="-log10 Poisson p-value")
    ax2.axhline(score_threshold, color="black", ls=":", label="Score threshold")
    ax2.set_ylabel("Enrichment score")
    ax2.set_xlabel("Genomic position (bp)")
    ax2.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig("figures/demo.png", dpi=120)
    plt.close(fig)

    print("Simulated {} true peaks, called {} peaks.".format(len(true_peaks), len(peaks)))
    print("Wrote figures/demo.png and results/summary.csv")


if __name__ == "__main__":
    main()
