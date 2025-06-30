import pytest
import pandas as pd
import gzip
from pathlib import Path
from plantvarfilter.parser import (
    smart_open,
    read_gene_traits
)

# ==== Tests for smart_open ====

def test_smart_open_reads_uncompressed(tmp_path):
    test_file = tmp_path / "test.csv"
    test_file.write_text("Gene_ID,Trait\nGene1,TraitA")
    with smart_open(test_file) as f:
        content = f.read()
    assert "Gene1" in content

def test_smart_open_reads_gz(tmp_path):
    test_file = tmp_path / "test.csv.gz"
    with gzip.open(test_file, "wt") as f:
        f.write("Gene_ID,Trait\nGene1,TraitA")
    with smart_open(test_file) as f:
        content = f.read()
    assert "Gene1" in content

def test_smart_open_raises_on_invalid_input():
    with pytest.raises(ValueError):
        smart_open(12345)  # Invalid input type

# ==== Tests for read_gene_traits ====

def test_read_gene_traits_handles_tsv_and_rename(tmp_path):
    tsv_file = tmp_path / "traits.tsv"
    tsv_file.write_text("Gene_ID\tTraitScore\nGene1\t0.9")
    df = read_gene_traits(tsv_file)
    assert "Gene" in df.columns
    assert df["Gene"].iloc[0] == "Gene1"
    assert df["TraitScore"].iloc[0] == 0.9

def test_read_gene_traits_handles_csv(tmp_path):
    csv_file = tmp_path / "traits.csv"
    csv_file.write_text("Gene,TraitScore\nGene1,0.9")
    df = read_gene_traits(csv_file)
    assert "Gene" in df.columns
    assert df["TraitScore"].iloc[0] == 0.9
