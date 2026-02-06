import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def _parse_marker(marker: str):
    """
    Extract chromosome and position from marker string.
    Examples:
      chr06:13755:C:T
      chr6:13755
      6:13755
    """
    if not isinstance(marker, str):
        return None, None

    m = re.match(r"(?:chr)?(\d+)[^\d]+(\d+)", marker)
    if not m:
        return None, None

    try:
        return int(m.group(1)), int(m.group(2))
    except Exception:
        return None, None


def plot_pangwas_manhattan(
    results_csv: str,
    out_png: str,
    p_col: str = "p_value",
    marker_col: str = "marker",
    alpha: float = 0.05,
    title: str = "panGWAS Manhattan plot",
):
    # ---------- basic checks ----------
    if not os.path.exists(results_csv):
        raise FileNotFoundError(results_csv)

    df = pd.read_csv(results_csv)

    if df.empty:
        print("[panGWAS][PLOT] CSV is empty, no Manhattan plot generated")
        return None

    if p_col not in df.columns or marker_col not in df.columns:
        raise ValueError(
            f"Required columns not found: {p_col}, {marker_col}"
        )

    # ---------- parse marker ----------
    parsed = df[marker_col].apply(_parse_marker)
    df["chrom"] = parsed.apply(lambda x: x[0])
    df["pos"] = parsed.apply(lambda x: x[1])

    df = df.dropna(subset=["chrom", "pos", p_col])

    if df.empty:
        print("[panGWAS][PLOT] No valid markers after parsing, skipping plot")
        return None

    # ---------- prepare data ----------
    df[p_col] = (
        pd.to_numeric(df[p_col], errors="coerce")
        .clip(lower=1e-300)
    )
    df = df.dropna(subset=[p_col])

    if df.empty:
        print("[panGWAS][PLOT] No valid p-values after cleanup")
        return None

    df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)

    # ---------- cumulative positions ----------
    chroms = df["chrom"].unique()
    chrom_offsets = {}
    xticks = []
    xtick_labels = []

    current_offset = 0
    for c in chroms:
        chrom_offsets[c] = current_offset
        chrom_len = df.loc[df["chrom"] == c, "pos"].max()

        xticks.append(current_offset + chrom_len / 2)
        xtick_labels.append(str(c))

        current_offset += chrom_len

    df["x"] = df.apply(
        lambda r: r["pos"] + chrom_offsets.get(r["chrom"], 0),
        axis=1,
    )
    df["mlogp"] = -np.log10(df[p_col])

    # ---------- plot ----------
    plt.figure(figsize=(14, 6))
    colors = ["#4C72B0", "#DD8452"]

    for i, c in enumerate(chroms):
        sub = df[df["chrom"] == c]
        plt.scatter(
            sub["x"],
            sub["mlogp"],
            s=10,
            color=colors[i % 2],
            alpha=0.8,
            linewidths=0,
        )

    if alpha and alpha > 0:
        sig = -np.log10(alpha)
        plt.axhline(sig, linestyle="--", color="red", linewidth=1)

    plt.xticks(xticks, xtick_labels)
    plt.xlabel("Chromosome")
    plt.ylabel(r"$-\log_{10}(P)$")
    plt.title(title)

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

    print(f"[panGWAS][PLOT] Manhattan plot saved: {out_png}")
    return out_png
