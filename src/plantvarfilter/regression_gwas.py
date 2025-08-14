# src/plantvarfilter/regression_gwas.py

import os
import math
import logging
from pathlib import Path
from typing import Optional, Tuple, List

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt

# ----------------------- utilities -----------------------

def fdr_bh(pvals: np.ndarray) -> np.ndarray:
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(1, n + 1))
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q, dtype=float)
    inv = np.empty_like(order)
    inv[order] = np.arange(n)
    out = np.clip(q, 0.0, 1.0)
    return out[inv]

def genomic_inflation_lambda(pvals: np.ndarray) -> float:
    # Chi-square with 1 df under null: median is ~0.4549364
    pv = np.asarray(pvals, dtype=float)
    pv = pv[np.isfinite(pv) & (pv > 0) & (pv <= 1)]
    if pv.size == 0:
        return float("nan")
    chisq = stats.chi2.isf(pv, df=1)
    lam = np.median(chisq) / stats.chi2.isf(0.5, df=1)
    return float(lam)

def save_qq_plot(pvals: np.ndarray, out_png: Path, title: str = "QQ Plot of GWAS P-values"):
    p = np.asarray(pvals, dtype=float)
    p = p[np.isfinite(p) & (p > 0) & (p <= 1)]
    if p.size == 0:
        logging.warning("No valid p-values for QQ plot.")
        return
    m = p.size
    exp = -np.log10((np.arange(1, m + 1)) / (m + 1))
    obs = -np.log10(np.sort(p))
    plt.figure(figsize=(8, 8))
    plt.scatter(exp, obs, s=10)
    maxv = max(exp.max(), obs.max())
    plt.plot([0, maxv], [0, maxv], linestyle="--")
    plt.xlabel("Expected -log10(P)")
    plt.ylabel("Observed -log10(P)")
    plt.title(title)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png)
    plt.close()

# ----------------------- core math -----------------------

def _residualize_against(Z: np.ndarray, v: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Regress vector v on columns of Z and return residuals and MSE of the fit.
    If Z has shape (n, k) and v has shape (n,), residuals r = v - Z*beta.
    """
    if Z is None or Z.size == 0:
        r = v.copy()
        mse = float(np.nan)
        return r, mse
    # QR for stability
    Q, R = np.linalg.qr(Z, mode="reduced")
    beta = np.linalg.solve(R, Q.T @ v)
    r = v - Z @ beta
    df = max(Z.shape[0] - Z.shape[1], 1)
    mse = float((r @ r) / df)
    return r, mse

def _fast_univariate_with_covariates(Xs: np.ndarray, y: np.ndarray, Z: Optional[np.ndarray]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    For each SNP column in Xs (n x m), test y ~ Z + x_j (linear).
    Returns beta, se, p, n per SNP.
    Approach: residualize y and each x_j w.r.t covariates Z, then simple OLS.
    """
    n, m = Xs.shape
    # residualize y
    ry, mse_y = _residualize_against(Z, y)
    # residualize all SNPs via projection on Q (from QR of Z)
    if Z is not None and Z.size:
        Q, R = np.linalg.qr(Z, mode="reduced")
        X_proj = Xs - Q @ (Q.T @ Xs)
    else:
        X_proj = Xs

    beta = np.full(m, np.nan, dtype=float)
    se   = np.full(m, np.nan, dtype=float)
    pval = np.full(m, 1.0, dtype=float)
    nobs = np.full(m, 0, dtype=int)

    ry_mean = np.nanmean(ry)

    for j in range(m):
        x = X_proj[:, j]
        mask = np.isfinite(x) & np.isfinite(ry)
        nj = int(mask.sum())
        if nj < 3:
            continue
        xv = x[mask]
        yv = ry[mask]
        # subtract means (since no intercept after residualization wrt Z)
        xv = xv - xv.mean()
        yv = yv - yv.mean()
        ssx = float(np.dot(xv, xv))
        if ssx <= 0:
            continue
        b = float(np.dot(xv, yv) / ssx)
        # residual variance
        resid = yv - b * xv
        df = nj - 1
        if df <= 0:
            continue
        s2 = float(np.dot(resid, resid) / df)
        seb = math.sqrt(s2 / ssx)
        t  = b / seb if seb > 0 else np.nan
        p  = 2 * stats.t.sf(abs(t), df) if np.isfinite(t) else 1.0
        beta[j] = b
        se[j]   = seb
        pval[j] = p
        nobs[j] = nj

    return beta, se, pval, nobs

# ----------------------- public API -----------------------

def run_snp_gwas_matrix(
    X: pd.DataFrame,        # dosage matrix (variants x samples)
    meta: pd.DataFrame,     # meta with index = X.index and columns CHROM, POS
    y: pd.Series,           # phenotype indexed by sample
    covariates: Optional[pd.DataFrame],
    out_dir: str
) -> pd.DataFrame:
    """
    Perform per-SNP linear regression: y ~ covariates + dosage.
    X: rows = SNP, cols = samples. y index must match X columns.
    covariates: DataFrame (rows = samples, columns = covariates), may be None.
    """
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    # align samples
    common = [s for s in X.columns if s in y.index]
    if covariates is not None:
        covariates = covariates.loc[common] if not covariates.empty else None
    X = X[common]
    y = y.loc[common].astype(float)

    # build Z (intercept + covariates)
    if covariates is not None and not covariates.empty:
        Z = covariates.astype(float).values
        Z = np.column_stack([np.ones(len(common)), Z])
    else:
        Z = np.ones((len(common), 1))

    # run
    beta, se, p, n = _fast_univariate_with_covariates(X.values.T, y.values, Z)

    res = pd.DataFrame({
        "SNP": X.index.values,
        "Beta": beta,
        "SE": se,
        "P": p,
        "N": n
    })

    # join meta (CHROM, POS) if available
    if meta is not None and not meta.empty:
        meta2 = meta.reset_index().rename(columns={"index": "SNP"})
        res = res.merge(meta2[["SNP", "CHROM", "POS"]], on="SNP", how="left")

    # FDR and lambdaGC + QQ plot
    valid_p = res["P"].values[np.isfinite(res["P"].values)]
    lam = genomic_inflation_lambda(valid_p)
    res["FDR_BH"] = fdr_bh(res["P"].fillna(1.0).values)

    res.to_csv(out_path / "gwas_snp_results.csv", index=False)

    # top hits
    if (res["P"].notna()).any():
        res.sort_values("P", ascending=True).head(100).to_csv(out_path / "gwas_top_hits.csv", index=False)

    # QQ plot
    save_qq_plot(valid_p, out_path / "plots" / "qq_plot.png")
    with open(out_path / "lambda_gc.txt", "w") as f:
        f.write(f"{lam:.4f}\n")
    logging.info(f"λGC = {lam:.4f} (saved to lambda_gc.txt)")

    return res

# Backward-compatible shim name
def run_regression_gwas_dynamic(*args, **kwargs):
    raise NotImplementedError("This module now provides SNP-level GWAS via run_snp_gwas_matrix().")

def run_regression_gwas(*args, **kwargs):
    raise NotImplementedError("Use run_snp_gwas_matrix() for SNP-level analysis.")
