import os
import math
from typing import Optional, Dict, List, Tuple, Any

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


PANGWAS_ALGORITHMS: Dict[str, Dict[str, Dict[str, str]]] = {
    "baseline": {
        "ttest": {
            "label": "Baseline • T-test (Welch)",
            "desc": "Two-sample Welch t-test between PAV=1 vs PAV=0."
        },
        "wilcoxon": {
            "label": "Baseline • Wilcoxon (Mann–Whitney U)",
            "desc": "Non-parametric rank test between PAV=1 vs PAV=0."
        },
        "fisher": {
            "label": "Baseline • Fisher Exact (Binary Trait)",
            "desc": "Exact test on 2x2 table (PAV vs binary phenotype)."
        },
    },
    "glm": {
        "ols": {
            "label": "GLM • Linear Regression (OLS) + covariates",
            "desc": "y ~ PAV + covariates (continuous trait)."
        },
        "logistic": {
            "label": "GLM • Logistic Regression + covariates",
            "desc": "logit(y) ~ PAV + covariates (binary trait)."
        },
    },
    "lmm": {
        "emmax": {
            "label": "LMM • Kinship Mixed Model (EMMAX)",
            "desc": "Linear mixed model with kinship correction (EMMAX-style)."
        }
    }
}


def list_pangwas_algorithms() -> List[Tuple[str, str, str]]:
    out: List[Tuple[str, str, str]] = []
    for cat, methods in PANGWAS_ALGORITHMS.items():
        for m, meta in methods.items():
            out.append((cat, m, meta["label"]))
    return out


def _read_table_auto(path: str) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)


def load_phenotype(pheno_path: str, trait: str) -> pd.Series:
    if not pheno_path or not os.path.exists(pheno_path):
        raise FileNotFoundError(f"Phenotype file not found: {pheno_path}")

    df = _read_table_auto(pheno_path)

    if "IID" in df.columns:
        sid = df["IID"].astype(str)
    else:
        sid = df.iloc[:, 0].astype(str)

    if trait not in df.columns:
        raise ValueError(f"Trait '{trait}' not found. Available: {list(df.columns)}")

    y = pd.to_numeric(df[trait], errors="coerce")
    return pd.Series(y.values, index=sid.values, name=trait).dropna()


def load_covariates(covar_path: str) -> pd.DataFrame:
    if not covar_path or not os.path.exists(covar_path):
        raise FileNotFoundError(f"Covariates file not found: {covar_path}")

    df = _read_table_auto(covar_path)

    if "IID" in df.columns:
        sid = df["IID"].astype(str)
        cov = df.drop(columns=["IID"])
    else:
        sid = df.iloc[:, 0].astype(str)
        cov = df.iloc[:, 1:]

    cov = cov.apply(pd.to_numeric, errors="coerce")
    cov.index = sid.values
    return cov.dropna(axis=0, how="any")


def _norm_sf(z: float) -> float:
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def _pval_from_z(z: float) -> float:
    z = abs(float(z))
    p = 2.0 * _norm_sf(z)
    return min(max(p, 0.0), 1.0)


def _rankdata(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a)
    order = np.argsort(a, kind="mergesort")
    ranks = np.empty(len(a), dtype=float)
    ranks[order] = np.arange(1, len(a) + 1, dtype=float)

    sorted_a = a[order]
    i = 0
    while i < len(sorted_a):
        j = i
        while j + 1 < len(sorted_a) and sorted_a[j + 1] == sorted_a[i]:
            j += 1
        if j > i:
            avg = 0.5 * (i + 1 + j + 1)
            ranks[order[i:j + 1]] = avg
        i = j + 1

    return ranks


def _align_inputs(
    pav: pd.DataFrame,
    y: pd.Series,
    covariates: Optional[pd.DataFrame] = None
):
    common = pav.index.astype(str).intersection(y.index.astype(str))
    if covariates is not None:
        common = common.intersection(covariates.index.astype(str))

    common = list(common)
    if len(common) < 5:
        raise ValueError(f"Not enough overlapping samples. Overlap={len(common)}")

    Xpav = pav.loc[common]
    yy = y.loc[common].astype(float).values
    C = covariates.loc[common].values.astype(float) if covariates is not None else None
    return Xpav, yy, C


def pangwas_ttest(pav: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
    print("[panGWAS][DEBUG] ENTER pangwas_ttest")

    print(f"[panGWAS][DEBUG] pav shape = {pav.shape}")
    print(f"[panGWAS][DEBUG] phenotype size = {len(y)}")

    X, yy, _ = _align_inputs(pav, y)

    print(f"[panGWAS][DEBUG] After align: X={X.shape}, y={yy.shape}")

    G = X.values.astype(np.int8)
    rows = []
    total = X.shape[1]

    print(f"[panGWAS] T-test started on {total} variants")

    for j, marker in enumerate(X.columns):
        if j % 10 == 0:
            print(f"[panGWAS][DEBUG] Variant {j}/{total}: {marker}")

        g = G[:, j]
        grp1 = yy[g == 1]
        grp0 = yy[g == 0]

        if len(grp1) < 2 or len(grp0) < 2:
            print(f"[panGWAS][DEBUG] Skip {marker} (n1={len(grp1)}, n0={len(grp0)})")
            continue

        stat, p = ttest_ind(grp1, grp0, equal_var=False, nan_policy="omit")

        rows.append((
            marker,
            float(p),
            float(grp1.mean() - grp0.mean()),
            len(grp1),
            len(grp0)
        ))

    print(f"[panGWAS] T-test finished, results={len(rows)} rows")

    return (
        pd.DataFrame(
            rows,
            columns=["marker", "p_value", "effect", "n1", "n0"]
        )
        .sort_values("p_value")
        .reset_index(drop=True)
    )


def pangwas_wilcoxon(pav: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
    X, yy, _ = _align_inputs(pav, y)

    rows = []
    for marker in X.columns:
        g = X[marker].values.astype(int)
        grp1, grp0 = yy[g == 1], yy[g == 0]
        if len(grp1) < 2 or len(grp0) < 2:
            continue

        allv = np.concatenate([grp1, grp0])
        ranks = _rankdata(allv)
        r1 = ranks[:len(grp1)].sum()
        U = r1 - len(grp1) * (len(grp1) + 1) / 2
        mu = len(grp1) * len(grp0) / 2
        sigma = math.sqrt(len(grp1) * len(grp0) * (len(grp1) + len(grp0) + 1) / 12)
        if sigma == 0:
            continue

        z = (U - mu) / sigma
        rows.append((
            marker,
            _pval_from_z(z),
            float(np.median(grp1) - np.median(grp0)),
            int(len(grp1)),
            int(len(grp0))
        ))

    return (
        pd.DataFrame(rows, columns=["marker", "p_value", "effect", "n1", "n0"])
        .sort_values("p_value")
        .reset_index(drop=True)
    )


def pangwas_fisher_exact(pav: pd.DataFrame, y: pd.Series) -> pd.DataFrame:
    X, yy, _ = _align_inputs(pav, y)
    uniq = sorted(set(yy))
    if len(uniq) != 2:
        raise ValueError("Fisher test requires binary phenotype.")

    yy01 = (yy == uniq[1]).astype(int)
    rows = []

    for marker in X.columns:
        g = X[marker].values.astype(int)
        a = int(np.sum((g == 1) & (yy01 == 1)))
        b = int(np.sum((g == 1) & (yy01 == 0)))
        c = int(np.sum((g == 0) & (yy01 == 1)))
        d = int(np.sum((g == 0) & (yy01 == 0)))
        if min(a + b, c + d) < 2:
            continue

        odds = ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
        rows.append((marker, float(odds)))

    df = pd.DataFrame(rows, columns=["marker", "odds_ratio"])
    df["p_value"] = np.nan
    return df


def pangwas_glm_ols(pav, y, covariates=None):
    X, yy, C = _align_inputs(pav, y, covariates)
    n = len(yy)
    C = np.zeros((n, 0)) if C is None else C

    rows = []
    for marker in X.columns:
        m = X[marker].values.reshape(-1, 1)
        D = np.column_stack([np.ones((n, 1)), m, C])
        beta, *_ = np.linalg.lstsq(D, yy, rcond=None)
        resid = yy - D @ beta
        s2 = (resid @ resid) / max(n - D.shape[1], 1)
        XtX_inv = np.linalg.pinv(D.T @ D)
        se = math.sqrt(max(0.0, s2 * XtX_inv[1, 1]))
        if se == 0:
            continue
        z = beta[1] / se
        rows.append((marker, _pval_from_z(z), float(beta[1]), float(se), int(n)))

    return (
        pd.DataFrame(rows, columns=["marker", "p_value", "beta", "se", "n"])
        .sort_values("p_value")
        .reset_index(drop=True)
    )


def pangwas_glm_logistic(pav, y, covariates=None):
    X, yy_raw, C = _align_inputs(pav, y, covariates)
    uniq = sorted(set(yy_raw))
    if len(uniq) != 2:
        raise ValueError("Logistic regression requires binary phenotype.")

    yy = (yy_raw == uniq[1]).astype(float)
    n = len(yy)
    C = np.zeros((n, 0)) if C is None else C

    def sigmoid(z):
        return 1.0 / (1.0 + np.exp(-np.clip(z, -30, 30)))

    rows = []
    for marker in X.columns:
        m = X[marker].values.reshape(-1, 1)
        D = np.column_stack([np.ones((n, 1)), m, C])
        beta = np.zeros(D.shape[1])

        for _ in range(50):
            mu = sigmoid(D @ beta)
            W = np.clip(mu * (1 - mu), 1e-8, None)
            z = D @ beta + (yy - mu) / W
            Dw = D * np.sqrt(W[:, None])
            zw = z * np.sqrt(W)
            beta_new, *_ = np.linalg.lstsq(Dw, zw, rcond=None)
            if np.max(np.abs(beta_new - beta)) < 1e-6:
                beta = beta_new
                break
            beta = beta_new

        XtWX = D.T @ (D * W[:, None])
        covb = np.linalg.pinv(XtWX)
        se = math.sqrt(max(0.0, covb[1, 1]))
        if se == 0:
            continue
        zstat = beta[1] / se
        rows.append((marker, _pval_from_z(zstat), float(beta[1]), float(se), int(n)))

    return (
        pd.DataFrame(rows, columns=["marker", "p_value", "beta", "se", "n"])
        .sort_values("p_value")
        .reset_index(drop=True)
    )


def compute_kinship_from_pav(pav):
    X = pav.values.astype(float)
    X -= X.mean(axis=0)
    std = X.std(axis=0, ddof=1)
    std[std == 0] = 1.0
    X /= std
    return (X @ X.T) / X.shape[1]


def pangwas_lmm_emmax(pav, y):
    X, yy, _ = _align_inputs(pav, y)
    K = compute_kinship_from_pav(X)
    w, U = np.linalg.eigh(K)
    w[w < 0] = 0
    delta = np.median(w)
    Vinv_sqrt = (U / np.sqrt(w + delta)) @ U.T
    yw = Vinv_sqrt @ yy
    X0 = Vinv_sqrt @ np.ones((len(yy), 1))

    rows = []
    for marker in X.columns:
        gw = Vinv_sqrt @ X[marker].values.reshape(-1, 1)
        D = np.column_stack([X0, gw])
        beta, *_ = np.linalg.lstsq(D, yw, rcond=None)
        resid = yw - D @ beta
        s2 = (resid @ resid) / max(len(yy) - 2, 1)
        XtX_inv = np.linalg.pinv(D.T @ D)
        se = math.sqrt(max(0.0, s2 * XtX_inv[1, 1]))
        if se == 0:
            continue
        z = beta[1] / se
        rows.append((marker, _pval_from_z(z), float(beta[1]), float(se)))

    return (
        pd.DataFrame(rows, columns=["marker", "p_value", "effect", "se"])
        .sort_values("p_value")
        .reset_index(drop=True)
    )


def run_pangwas(
    pav,
    y,
    category,
    method,
    covariates=None,
    kinship=None,
    extra=None
):
    print("[panGWAS][DEBUG] ENTER run_pangwas")

    category = str(category).lower().strip()
    method = str(method).lower().strip()

    print(f"[panGWAS][DEBUG] category={category}, method={method}")
    print(f"[panGWAS][DEBUG] pav shape={getattr(pav, 'shape', None)}")
    print(f"[panGWAS][DEBUG] phenotype size={len(y) if y is not None else None}")

    try:
        if category == "baseline":
            if method == "ttest":
                return pangwas_ttest(pav, y)
            if method == "wilcoxon":
                return pangwas_wilcoxon(pav, y)
            if method == "fisher":
                return pangwas_fisher_exact(pav, y)

        if category == "glm":
            if method == "ols":
                return pangwas_glm_ols(pav, y, covariates)
            if method == "logistic":
                return pangwas_glm_logistic(pav, y, covariates)

        if category == "lmm":
            if method == "emmax":
                return pangwas_lmm_emmax(pav, y)

        raise ValueError(f"Unknown panGWAS configuration: {category}/{method}")

    except Exception as e:
        print("[panGWAS][ERROR] run_pangwas failed:", repr(e))
        raise
