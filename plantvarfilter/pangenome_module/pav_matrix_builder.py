import os
from typing import List, Optional, Set, Tuple

import numpy as np
import pandas as pd


class GenePresenceAbsenceBuilder:
    def __init__(
        self,
        gff_dir: str,
        sample_list: Optional[str] = None,
        gene_id_field: str = "ID",
        max_genes: Optional[int] = None,
        sort_genes: bool = True,
    ):
        self.gff_dir = gff_dir
        self.sample_list = sample_list
        self.gene_id_field = gene_id_field
        self.max_genes = max_genes
        self.sort_genes = sort_genes

        self.samples = self._load_samples()

    def _load_samples(self) -> List[str]:
        if self.sample_list:
            with open(self.sample_list, "r", encoding="utf-8") as f:
                return [line.strip() for line in f if line.strip()]

        samples = []
        for fname in os.listdir(self.gff_dir):
            if fname.endswith(".gff"):
                samples.append(fname[:-4])
            elif fname.endswith(".gtf"):
                samples.append(fname[:-4])
        return samples

    def _find_gff(self, sample: str) -> str:
        for ext in (".gff", ".gtf"):
            path = os.path.join(self.gff_dir, sample + ext)
            if os.path.exists(path):
                return path
        raise FileNotFoundError(f"No GFF/GTF found for sample: {sample}")

    def _extract_gene_id(self, attributes: str) -> Optional[str]:
        key = self.gene_id_field + "="
        for field in attributes.split(";"):
            field = field.strip()
            if field.startswith(key):
                val = field.split("=", 1)[1].strip()
                return val if val else None
        return None

    def _parse_gff(self, gff_path: str) -> Set[str]:
        genes: Set[str] = set()
        with open(gff_path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line or line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 9:
                    continue
                if parts[2] != "gene":
                    continue

                gid = self._extract_gene_id(parts[8])
                if gid:
                    genes.add(gid)
                    if self.max_genes is not None and len(genes) >= self.max_genes:
                        break
        return genes

    def build_matrix(self) -> pd.DataFrame:
        records: List[Tuple[str, str]] = []
        gene_seen: Set[str] = set()

        for sample in self.samples:
            gff_path = self._find_gff(sample)
            genes = self._parse_gff(gff_path)

            for gene in genes:
                records.append((sample, gene))
                gene_seen.add(gene)

            if self.max_genes is not None and len(gene_seen) >= self.max_genes:
                break

        if not records:
            out = pd.DataFrame(index=pd.Index(self.samples, name="Sample"))
            return out

        df = pd.DataFrame(records, columns=["Sample", "Gene"])
        df["PAV"] = 1

        pav = (
            df.pivot_table(
                index="Sample",
                columns="Gene",
                values="PAV",
                aggfunc="max",
                fill_value=0,
            )
            .astype(np.int8)
        )

        pav.index.name = "Sample"
        if self.sort_genes:
            pav = pav.sort_index(axis=1)

        return pav

    def save_matrix(self, out_path: str) -> str:
        df = self.build_matrix()
        df.to_csv(out_path)
        return out_path
