import pandas as pd
import pytest
from io import StringIO
from plantvarfilter.annotator import annotate_variants_with_genes, annotate_with_traits
from plantvarfilter.filter import improved_filter_variants


# ==== Tests for improved_filter_variants ====

def test_improved_filter_variants_reads_and_filters(tmp_path):
    vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1000\t.\tA\tT\t.\t.\tCSQ=T|missense_variant
1\t2000\t.\tG\tC\t.\t.\tCSQ=C|intron_variant
"""
    vcf_path = tmp_path / "test.vcf"
    with open(vcf_path, "w") as f:
        f.write(vcf_content)

    with open(vcf_path, "rb") as stream:
        output_path = improved_filter_variants(
            stream,
            include_intergenic=False,
            store_as_feather=True,
            consequence_types=["missense_variant"]
        )

    assert output_path.endswith(".feather")
    df = pd.read_feather(output_path)
    assert len(df) == 1
    assert df["POS"].iloc[0] == 1000
    assert "missense_variant" in df["Consequence"].iloc[0]

