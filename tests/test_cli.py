import pandas as pd
import pytest
from io import StringIO
from plantvarfilter.annotator import annotate_variants_with_genes, annotate_with_traits
from plantvarfilter.filter import improved_filter_variants
from plantvarfilter.parser import parse_info_field, parse_csq_field, smart_open, read_gene_traits
import subprocess
import sys
from pathlib import Path

def test_cli_runs_successfully(tmp_path):
    vcf_path = tmp_path / "test.vcf"
    gff_path = tmp_path / "test.gff"
    traits_path = tmp_path / "traits.csv"
    config_path = tmp_path / "config.json"
    output_dir = tmp_path / "cli_output"
    output_dir.mkdir()

    vcf_path.write_text("""##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t1000\t.\tA\tT\t.\t.\tCSQ=T|missense_variant
""")

    gff_path.write_text("""##gff-version 3
1\tsource\tgene\t900\t1100\t.\t+\t.\tID=Gene1
""")

    traits_path.write_text("Gene,Trait_Score\nGene1,0.85")

    config = {
        "vcf": str(vcf_path),
        "gff": str(gff_path),
        "traits": str(traits_path),
        "include_intergenic": True,
        "consequence_types": ["missense_variant"],
        "output_format": "csv",
        "output": str(output_dir / "filtered_variants.csv"),
        "plot": False,
        "gwas": False,
        "output_dir": str(output_dir)
    }

    import json
    config_path.write_text(json.dumps(config, indent=2))

    result = subprocess.run(
        [sys.executable, "-m", "plantvarfilter", "run", "--config", str(config_path)],
        capture_output=True,
        text=True
    )

    print(result.stdout)
    print(result.stderr)

    assert result.returncode == 0
    assert (output_dir / "filtered_variants.csv").exists()


