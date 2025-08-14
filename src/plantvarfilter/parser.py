import pandas as pd
import gzip
from typing import Union, TextIO
from pathlib import Path
import re
from io import StringIO

MISSING_CODES = {"-9", "NA", "NaN", ".", ""}


def smart_open(file_or_path: Union[str, Path, TextIO]) -> TextIO:
    """Open text file (optionally .gz) and return a text stream."""
    if isinstance(file_or_path, (str, Path)):
        p = str(file_or_path)
        if p.endswith(".gz"):
            return gzip.open(p, "rt", encoding="utf-8")
        return open(p, "r", encoding="utf-8")
    if hasattr(file_or_path, "read"):
        return file_or_path
    raise ValueError("Unsupported file input type")


def _parse_single_col_traits(lines: list[str]) -> pd.DataFrame:
    """
    Parse free-text lines like: 'V002 V002 73.18' or 'V002,73.18'.
    Returns a DataFrame with columns ['Sample', 'Trait_Score'].
    """
    samples, scores = [], []
    for raw in lines:
        s = str(raw).strip().replace(",", " ")
        if not s:
            continue
        toks = s.split()
        if len(toks) < 2:
            continue
        sample_token = next((t for t in toks if re.fullmatch(r"[Vv]\d+", t)), toks[0])
        val = None
        for t in reversed(toks):
            try:
                val = float(t)
                break
            except Exception:
                continue
        if val is None or str(val) in MISSING_CODES or val == -9:
            continue
        samples.append(str(sample_token).upper().strip())
        scores.append(val)
    if not samples:
        return pd.DataFrame(columns=["Sample", "Trait_Score"])
    df = pd.DataFrame({"Sample": samples, "Trait_Score": scores})
    df["Sample"] = df["Sample"].astype(str).str.strip().str.upper()
    return df.dropna(subset=["Trait_Score"])


def _read_csv_smart_text(content: str) -> pd.DataFrame:
    """
    Try reading CSV/TSV with automatic delimiter inference, then tab, then whitespace.
    Returns a pandas DataFrame (no column normalization here).
    """
    try:
        return pd.read_csv(StringIO(content), sep=None, engine="python", dtype=str)
    except Exception:
        pass
    try:
        return pd.read_csv(StringIO(content), sep="\t", dtype=str)
    except Exception:
        pass
    return pd.read_csv(StringIO(content), sep=r"\s+", engine="python", dtype=str, header=None)


def _normalize_traits_df(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize a raw traits DataFrame into either:
      ['Sample','Trait_Score'] or ['Gene','Trait_Score'].
    Handles:
      - duplicated headers
      - header mis-detection
      - alias columns
      - triple form: name/name/value
      - single text column fallback
    """
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]

    lower_cols = [str(c).lower() for c in df.columns]
    looks_like_data_header = False
    if len(df) > 0:
        row0 = df.iloc[0].fillna("").tolist()
        looks_like_data_header = all(
            bool(re.fullmatch(r"[Vv]\d+", str(x))) or str(x).replace(".", "", 1).isdigit()
            for x in row0
        )
    if not any(c in lower_cols for c in ("sample", "gene", "id", "name")) and looks_like_data_header:
        df = pd.read_csv(StringIO(df.to_csv(index=False, header=False)), header=None)

    cols_lower = {c.lower(): c for c in df.columns.astype(str)}

    if "sample" not in cols_lower and "gene" not in cols_lower:
        for a in ("id", "name"):
            if a in cols_lower:
                df.rename(columns={cols_lower[a]: "Sample"}, inplace=True)
                cols_lower["sample"] = "Sample"
                break

    if "trait_score" not in cols_lower:
        for a in ("score", "phenotype", "value", "trait", "pheno", "y"):
            if a in cols_lower:
                df.rename(columns={cols_lower[a]: "Trait_Score"}, inplace=True)
                cols_lower["trait_score"] = "Trait_Score"
                break

    if df.shape[1] >= 3:
        c0, c1 = df.columns[0], df.columns[1]
        try:
            same = (
                df[c0].astype(str).str.strip().str.upper()
                == df[c1].astype(str).str.strip().str.upper()
            )
            if same.all():
                out = pd.DataFrame(
                    {
                        "Sample": df.iloc[:, 0].astype(str).str.strip().str.upper(),
                        "Trait_Score": pd.to_numeric(df.iloc[:, 2], errors="coerce"),
                    }
                )
                out["Trait_Score"] = out["Trait_Score"].where(~out["Trait_Score"].isin({-9}), pd.NA)
                out = out.replace({"Trait_Score": {m: pd.NA for m in MISSING_CODES}})
                return out.dropna(subset=["Trait_Score"])
        except Exception:
            pass

    if df.shape[1] == 1:
        return _parse_single_col_traits(df.iloc[:, 0].astype(str).tolist())

    df = df.replace({m: pd.NA for m in MISSING_CODES})

    if "Sample" in df.columns:
        df["Sample"] = df["Sample"].astype(str).str.strip().str.upper()
        target = None
        if "Trait_Score" in df.columns:
            target = "Trait_Score"
        else:
            numeric_candidates = []
            for c in df.columns:
                if c == "Sample":
                    continue
                col_num = pd.to_numeric(df[c], errors="coerce")
                if col_num.notna().sum() > 0:
                    numeric_candidates.append(c)
            if numeric_candidates:
                target = numeric_candidates[0]
                if target != "Trait_Score":
                    df.rename(columns={target: "Trait_Score"}, inplace=True)
                    target = "Trait_Score"
        if target is None:
            return _parse_single_col_traits(df.astype(str).agg(" ".join, axis=1).tolist())
        df["Trait_Score"] = pd.to_numeric(df["Trait_Score"], errors="coerce")
        return df[["Sample", "Trait_Score"]].dropna(subset=["Trait_Score"])

    if "Gene" in df.columns:
        df["Gene"] = df["Gene"].astype(str).str.strip()
        if "Trait_Score" not in df.columns:
            numeric_candidates = []
            for c in df.columns:
                if c == "Gene":
                    continue
                col_num = pd.to_numeric(df[c], errors="coerce")
                if col_num.notna().sum() > 0:
                    numeric_candidates.append(c)
            if numeric_candidates:
                df.rename(columns={numeric_candidates[0]: "Trait_Score"}, inplace=True)
        df["Trait_Score"] = pd.to_numeric(df["Trait_Score"], errors="coerce")
        return df[["Gene", "Trait_Score"]].dropna(subset=["Trait_Score"])

    return _parse_single_col_traits(df.astype(str).agg(" ".join, axis=1).tolist())


def read_gene_traits(file_or_path: Union[str, TextIO]) -> pd.DataFrame:
    """
    Smart traits reader that always returns either:
      ['Sample', 'Trait_Score']  or  ['Gene', 'Trait_Score'].
    Supports:
      - CSV/TSV with alias headers and automatic delimiter inference
      - free-text lines like 'V002 V002 73.18'
      - header mis-detection fallback
    Missing codes (-9, NA, NaN, .) are treated as NA.
    """
    try:
        with smart_open(file_or_path) as fh:
            content = fh.read()
    except Exception as e:
        raise RuntimeError(f"Failed to open traits file: {e}")

    if "\n" in content and ("," not in content and "\t" not in content):
        lines = [ln for ln in content.splitlines() if ln.strip()]
        hits = sum(bool(re.search(r"[Vv]\d+", ln)) for ln in lines)
        if hits >= max(3, int(0.6 * len(lines))):
            return _parse_single_col_traits(lines)

    df = _read_csv_smart_text(content)
    return _normalize_traits_df(df)
