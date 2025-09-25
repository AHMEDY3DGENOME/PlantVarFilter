# annotation_utils.py
import os
import gzip
import io
from typing import Dict, Tuple, Optional

import numpy as np
import pandas as pd


class Annotator:
    """
    Lightweight GWAS annotator:
      - Loads GTF/GFF (optionally .gz)
      - Extracts 'gene' features with gene_id / gene_name
      - Builds per-chromosome sorted arrays for fast nearest-gene lookup
      - Annotates GWAS table with nearest gene within a user window (kb)
    """

    def __init__(self):
        self.genes_df: Optional[pd.DataFrame] = None
        self.index: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}

    # ---------- I/O helpers ----------
    @staticmethod
    def _open_any(path: str):
        if str(path).endswith(".gz"):
            with gzip.open(path, "rb") as fh:
                data = fh.read()
            return io.StringIO(data.decode("utf-8", errors="replace"))
        return open(path, "r", encoding="utf-8", errors="replace")

    @staticmethod
    def _parse_attrs(attr: str) -> Dict[str, str]:
        out = {}
        if not isinstance(attr, str):
            return out

        # GTF style: key "value"; key "value";
        # GFF style: key=value;key2=value2
        parts = []
        # try GFF k=v
        if "=" in attr and ";" in attr:
            parts = [p for p in attr.strip().split(";") if p.strip()]
            for p in parts:
                if "=" in p:
                    k, v = p.split("=", 1)
                    out[k.strip()] = v.strip().strip('"')
        else:
            # fall back to GTF
            semi = [p for p in attr.strip().split(";") if p.strip()]
            for p in semi:
                tokens = p.strip().split()
                if len(tokens) >= 2:
                    k = tokens[0].strip()
                    v = " ".join(tokens[1:]).strip().strip('"')
                    out[k] = v
        return out

    # ---------- GTF/GFF load ----------
    def load_gtf_or_gff(self, path: str) -> pd.DataFrame:
        cols = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
        rows = []
        with self._open_any(path) as fh:
            for ln in fh:
                if not ln or ln.startswith("#"):
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                seq, src, feat, st, en, sc, strand, frame, attr = parts[:9]
                if feat.lower() != "gene":
                    continue
                try:
                    st_i = int(st)
                    en_i = int(en)
                except Exception:
                    continue
                a = self._parse_attrs(attr)
                gene_id = a.get("gene_id") or a.get("ID") or a.get("GeneID") or ""
                gene_name = a.get("gene_name") or a.get("Name") or a.get("gene") or gene_id
                rows.append((seq, st_i, en_i, gene_id, gene_name, strand))

        if not rows:
            raise ValueError("No 'gene' features found in supplied GTF/GFF.")

        df = pd.DataFrame(rows, columns=["Chr", "Start", "End", "GeneID", "GeneName", "Strand"])
        # Normalize chromosome naming (remove 'chr' prefix if mixed)
        if df["Chr"].str.startswith("chr").any():
            df["Chr"] = df["Chr"].str.replace("^chr", "", regex=True)
        self.genes_df = df
        return df

    # ---------- Index build ----------
    def build_index(self):
        if self.genes_df is None or self.genes_df.empty:
            raise RuntimeError("Call load_gtf_or_gff first.")

        self.index.clear()
        for chrom, sub in self.genes_df.groupby("Chr", sort=False):
            # sort by Start for binary search
            sub = sub.sort_values("Start").reset_index(drop=True)
            starts = sub["Start"].to_numpy(dtype=np.int64)
            ends = sub["End"].to_numpy(dtype=np.int64)
            names = sub["GeneName"].astype(str).to_numpy()
            ids = sub["GeneID"].astype(str).to_numpy()
            self.index[str(chrom)] = (starts, ends, names, ids)

    # ---------- nearest gene ----------
    @staticmethod
    def _nearest_distance(pos: int, starts: np.ndarray, ends: np.ndarray) -> Tuple[int, int]:
        """
        Returns (idx, distance_bp). If pos is inside a gene interval, distance=0.
        Otherwise min distance to closest interval edge. idx is index into arrays.
        """
        if starts.size == 0:
            return -1, np.iinfo(np.int64).max

        # locate insertion point
        i = np.searchsorted(starts, pos, side="left")

        candidates = []
        if i < starts.size:
            # interval to the right
            d = 0 if (pos >= starts[i] and pos <= ends[i]) else (starts[i] - pos if pos < starts[i] else pos - ends[i])
            candidates.append((i, abs(int(d))))
        if i > 0:
            j = i - 1
            d = 0 if (pos >= starts[j] and pos <= ends[j]) else (starts[j] - pos if pos < starts[j] else pos - ends[j])
            candidates.append((j, abs(int(d))))

        idx, dist = min(candidates, key=lambda t: t[1]) if candidates else (-1, np.iinfo(np.int64).max)
        return idx, dist

    # ---------- public API ----------
    def annotate_gwas(
        self,
        gwas_df: pd.DataFrame,
        window_kb: int = 50,
        chr_col: str = "Chr",
        pos_col: str = "ChrPos",
    ) -> pd.DataFrame:
        """
        Returns a copy of gwas_df with columns:
          - NearestGene
          - NearestGeneID
          - DistanceToGene
          - InWindow  (True if distance <= window_kb*1000)
        """
        if self.genes_df is None or not self.index:
            raise RuntimeError("Load GTF/GFF and call build_index() before annotate_gwas().")

        df = gwas_df.copy()
        if df[chr_col].dtype != object:
            df[chr_col] = df[chr_col].astype(str)
        # normalize 'chr' prefix as index uses normalized strings
        df[chr_col] = df[chr_col].str.replace("^chr", "", regex=True)

        nearest_gene = []
        nearest_id = []
        distance_bp = []
        in_window = []
        max_bp = int(window_kb) * 1000

        for chrom, pos in zip(df[chr_col].tolist(), df[pos_col].tolist()):
            tup = self.index.get(str(chrom))
            if tup is None:
                nearest_gene.append("")
                nearest_id.append("")
                distance_bp.append(np.iinfo(np.int64).max)
                in_window.append(False)
                continue
            starts, ends, names, ids = tup
            idx, dist = self._nearest_distance(int(pos), starts, ends)
            if idx >= 0:
                nearest_gene.append(str(names[idx]))
                nearest_id.append(str(ids[idx]))
                distance_bp.append(int(dist))
                in_window.append(dist <= max_bp)
            else:
                nearest_gene.append("")
                nearest_id.append("")
                distance_bp.append(np.iinfo(np.int64).max)
                in_window.append(False)

        df["NearestGene"] = nearest_gene
        df["NearestGeneID"] = nearest_id
        df["DistanceToGene"] = distance_bp
        df["InWindow"] = in_window
        return df

    def annotate_and_save(
        self,
        gwas_df: pd.DataFrame,
        out_csv: str,
        window_kb: int = 50,
        chr_col: str = "Chr",
        pos_col: str = "ChrPos",
    ) -> str:
        df_annot = self.annotate_gwas(gwas_df, window_kb=window_kb, chr_col=chr_col, pos_col=pos_col)
        out = os.path.splitext(out_csv)[0] + "_annotated.csv"
        df_annot.to_csv(out, index=False)
        return out

    # ---------- optional: region filter ----------
    @staticmethod
    def filter_region(df: pd.DataFrame, region: str, chr_col="Chr", pos_col="ChrPos") -> pd.DataFrame:
        """
        region format: 'chr:start-end' or 'chrom:start-end' (chr prefix optional)
        """
        if not region:
            return df
        region = region.strip()
        if ":" not in region or "-" not in region:
            return df
        chrom, rest = region.split(":", 1)
        start, end = rest.split("-", 1)
        chrom = chrom.replace("chr", "")
        start_i = int(start.replace(",", ""))
        end_i = int(end.replace(",", ""))
        df2 = df.copy()
        if df2[chr_col].dtype != object:
            df2[chr_col] = df2[chr_col].astype(str)
        df2[chr_col] = df2[chr_col].str.replace("^chr", "", regex=True)
        return df2[(df2[chr_col] == chrom) & (df2[pos_col].between(start_i, end_i))]
