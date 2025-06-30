import pandas as pd
import pytest
from plantvarfilter.annotator import annotate_variants_with_genes, annotate_with_traits

# ==== Helper Fixtures ====
@pytest.fixture
def gene_db():
    return {
        "Gene1": {"seqid": "1", "start": 900, "end": 1100, "strand": "+"},
        "Gene2": {"seqid": "1", "start": 2000, "end": 2200, "strand": "+"},
        "Gene3": {"seqid": "2", "start": 500, "end": 800, "strand": "-"},
    }

# ==== Tests for annotate_variants_with_genes ====

def test_match_within_gene(gene_db):
    vcf_df = pd.DataFrame({
        "CHROM": ["1"],
        "POS": [1000],
        "CSQ": ["T|missense_variant"],
        "Gene": [None]
    })
    result = annotate_variants_with_genes(vcf_df, gene_db)
    assert result["Gene"].iloc[0] == "Gene1"

def test_intergenic_within_same_chrom(gene_db):
    vcf_df = pd.DataFrame({
        "CHROM": ["1"],
        "POS": [3000],
        "CSQ": ["A|synonymous_variant"],
        "Gene": [None]
    })
    result = annotate_variants_with_genes(vcf_df, gene_db)
    assert result["Gene"].iloc[0] == "Gene2"  # Closest

def test_intergenic_excluded(gene_db):
    vcf_df = pd.DataFrame({
        "CHROM": ["1"],
        "POS": [3000],
        "CSQ": ["A|synonymous_variant"],
        "Gene": [None]
    })
    result = annotate_variants_with_genes(vcf_df, gene_db, include_intergenic=False)
    assert pd.isna(result["Gene"].iloc[0])

def test_different_chrom_returns_none(gene_db):
    vcf_df = pd.DataFrame({
        "CHROM": ["3"],
        "POS": [100],
        "CSQ": ["G|stop_gained"],
        "Gene": [None]
    })
    result = annotate_variants_with_genes(vcf_df, gene_db)
    assert pd.isna(result["Gene"].iloc[0])

# ==== Tests for annotate_with_traits ====

def test_trait_annotation_found():
    variants_df = pd.DataFrame({"Gene": ["Gene1"]})
    traits_df = pd.DataFrame({"Gene": ["Gene1"], "TraitScore": [0.85]})
    result = annotate_with_traits(variants_df, traits_df)
    assert result["TraitScore"].iloc[0] == 0.85

def test_trait_annotation_not_found():
    variants_df = pd.DataFrame({"Gene": ["GeneX"]})
    traits_df = pd.DataFrame({"Gene": ["Gene1"], "TraitScore": [0.85]})
    result = annotate_with_traits(variants_df, traits_df)
    assert pd.isna(result["TraitScore"].iloc[0])

def test_trait_annotation_keep_unmatched_false():
    variants_df = pd.DataFrame({"Gene": ["GeneX"]})
    traits_df = pd.DataFrame({"Gene": ["Gene1"], "TraitScore": [0.85]})
    result = annotate_with_traits(variants_df, traits_df, keep_unmatched=False)
    assert result.empty
