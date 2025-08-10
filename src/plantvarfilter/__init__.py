import os, json
from pathlib import Path

USER_DATA_DIR = Path.home() / ".plantvarfilter_data"
CONFIG_PATH = USER_DATA_DIR / "config.json"
EXAMPLE_CONFIG = {
    "vcf": "input/your_data.vcf.gz",
    "gff": "input/your_genes.gff3.gz",
    "traits": "input/your_traits.csv",
    "include_intergenic": True,
    "consequence_types": ["missense_variant", "stop_gained"],
    "output_format": "csv",
    "output": "filtered_output.csv",
    "plot": True,
    "gwas": True,
}

def initialize_user_data():
    if not USER_DATA_DIR.exists():
        (USER_DATA_DIR / "input").mkdir(parents=True, exist_ok=True)
        (USER_DATA_DIR / "output").mkdir(parents=True, exist_ok=True)
        with open(CONFIG_PATH, "w") as f:
            json.dump(EXAMPLE_CONFIG, f, indent=2)

from .regression_gwas import run_regression_gwas, run_regression_gwas_dynamic

__all__ = ["run_regression_gwas", "run_regression_gwas_dynamic", "initialize_user_data"]
