# A from-scratch sliding-window peak caller with Poisson background modelling.
import argparse
import numpy as np
from scipy.stats import poisson


def sliding_counts(coverage, window):
    # Sum coverage over a sliding window using a cumulative-sum trick.
    csum = np.concatenate(([0.0], np.cumsum(coverage.astype(float))))
    n = len(coverage)
    counts = np.zeros(n)
    half = window // 2
    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        counts[i] = csum[end] - csum[start]
    return counts


def local_background(coverage, flank):
    # Estimate a local mean rate per base using a wide flanking window.
    csum = np.concatenate(([0.0], np.cumsum(coverage.astype(float))))
    n = len(coverage)
    bg = np.zeros(n)
    half = flank // 2
    for i in range(n):
        start = max(0, i - half)
        end = min(n, i + half + 1)
        width = end - start
        bg[i] = (csum[end] - csum[start]) / max(width, 1)
    return bg


def compute_scores(coverage, window=50, flank=1000):
    # Return a -log10 Poisson survival score for every position.
    counts = sliding_counts(coverage, window)
    bg_rate = local_background(coverage, flank)
    global_rate = coverage.mean()
    lam = np.maximum(bg_rate, global_rate) * window
    lam = np.maximum(lam, 1e-9)
    # Probability of observing at least this many reads under background.
    pval = poisson.sf(counts - 1, lam)
    pval = np.clip(pval, 1e-300, 1.0)
    scores = -np.log10(pval)
    return scores


def merge_intervals(positions, coverage, scores, min_gap):
    # Merge nearby flagged positions into peaks and record summits.
    if len(positions) == 0:
        return []
    peaks = []
    start = positions[0]
    prev = positions[0]
    for p in positions[1:]:
        if p - prev <= min_gap:
            prev = p
        else:
            peaks.append((start, prev))
            start = p
            prev = p
    peaks.append((start, prev))
    results = []
    for s, e in peaks:
        region = coverage[s:e + 1]
        summit = s + int(np.argmax(region))
        peak_score = float(scores[s:e + 1].max())
        results.append({
            "start": int(s),
            "end": int(e),
            "summit": int(summit),
            "score": round(peak_score, 3),
            "max_coverage": float(coverage[summit]),
        })
    return results


def call_peaks(coverage, window=50, flank=1000, score_threshold=5.0, min_gap=100):
    # Main entry point: coverage is a 1D numpy array.
    coverage = np.asarray(coverage, dtype=float)
    scores = compute_scores(coverage, window=window, flank=flank)
    flagged = np.where(scores >= score_threshold)[0]
    peaks = merge_intervals(flagged, coverage, scores, min_gap=min_gap)
    return peaks, scores


def load_coverage(path):
    # Load a one-value-per-line coverage file.
    return np.loadtxt(path)


def main():
    parser = argparse.ArgumentParser(description="Sliding-window Poisson peak caller.")
    parser.add_argument("coverage", help="Path to one-value-per-line coverage file.")
    parser.add_argument("--window", type=int, default=50)
    parser.add_argument("--flank", type=int, default=1000)
    parser.add_argument("--score-threshold", type=float, default=5.0)
    parser.add_argument("--min-gap", type=int, default=100)
    parser.add_argument("--out", default="peaks.tsv")
    args = parser.parse_args()

    coverage = load_coverage(args.coverage)
    peaks, _ = call_peaks(
        coverage,
        window=args.window,
        flank=args.flank,
        score_threshold=args.score_threshold,
        min_gap=args.min_gap,
    )
    with open(args.out, "w") as fh:
        fh.write("start\tend\tsummit\tscore\tmax_coverage\n")
        for pk in peaks:
            fh.write("{start}\t{end}\t{summit}\t{score}\t{max_coverage}\n".format(**pk))
    print("Called {} peaks. Written to {}".format(len(peaks), args.out))


if __name__ == "__main__":
    main()
