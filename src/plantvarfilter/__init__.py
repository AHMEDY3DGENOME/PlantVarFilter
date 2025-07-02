import os
from pathlib import Path
import json

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
    "gwas": True
}

def initialize_user_data():
    """
    Create ~/.plantvarfilter_data/ with input/, output/, and default config.json if not exists.
    """
    if not USER_DATA_DIR.exists():
        os.makedirs(USER_DATA_DIR / "input", exist_ok=True)
        os.makedirs(USER_DATA_DIR / "output", exist_ok=True)
        with open(CONFIG_PATH, "w") as f:
            json.dump(EXAMPLE_CONFIG, f, indent=2)

#Import regression GWAS so it's available directly from the package
from .regression_gwas import run_regression_gwas
