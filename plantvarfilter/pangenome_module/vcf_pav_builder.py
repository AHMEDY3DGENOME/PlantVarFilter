import os
import numpy as np
import pandas as pd
from cyvcf2 import VCF


class VariantPAVBuilder:
    """
    Build Presence/Absence matrix (0/1) from a multi-sample VCF/VCF.GZ.
    Presence = 1 if genotype is not homozygous reference (0/0).
    Output shape: samples x variants
    """

    def __init__(self, vcf_path: str, max_variants: int = 0):
        self.vcf_path = vcf_path
        self.max_variants = int(max_variants or 0)

    def build_matrix(self) -> pd.DataFrame:
        if not self.vcf_path or not os.path.exists(self.vcf_path):
            raise FileNotFoundError(f"VCF not found: {self.vcf_path}")

        vcf = VCF(self.vcf_path)
        samples = vcf.samples

        if not samples:
            raise RuntimeError("No samples found in VCF")

        variant_ids = []
        rows = []

        for i, var in enumerate(vcf):
            if self.max_variants and i >= self.max_variants:
                break

            # genotypes: [ [a, b, phased], ... ]
            gts = var.genotypes
            row = np.array(
                [1 if (g[0] != 0 or g[1] != 0) else 0 for g in gts],
                dtype=np.int8,
            )

            vid = var.ID
            if not vid or vid == ".":
                vid = f"{var.CHROM}:{var.POS}:{var.REF}:{','.join(var.ALT)}"

            variant_ids.append(vid)
            rows.append(row)

        if not rows:
            return pd.DataFrame(index=pd.Index(samples, name="Sample"))

        data = np.vstack(rows).T
        df = pd.DataFrame(data, index=samples, columns=variant_ids)
        df.index.name = "Sample"
        return df
