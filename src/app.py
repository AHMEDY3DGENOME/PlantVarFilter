import streamlit as st
import time
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
import io
import pandas as pd

from filter import (
    improved_filter_variants as filter_variants,
    build_gene_db,
    annotate_variants_with_genes,
    annotate_with_traits
)
from parser import smart_open, read_gene_traits

st.set_page_config(page_title="PlantVarFilter", layout="wide")
st.title("🌿 PlantVarFilter: Variant Filtering for Plant Genomics")

method = st.radio("Select file input method:", ["Upload Files", "Use Local Paths"])

vcf_file = gff_file = traits_file = None
vcf_path = gff_path = traits_path = ""

include_intergenic = st.checkbox("Include intergenic variants", value=True)

# File upload or path input
if method == "Upload Files":
    vcf_file = st.file_uploader("Upload VCF/VCF.GZ File", type=["vcf", "vcf.gz"])
    gff_file = st.file_uploader("Upload GFF3/GFF3.GZ File", type=["gff", "gff3", "gff.gz", "gff3.gz"])
    traits_file = st.file_uploader("Upload Trait File (CSV/TSV)", type=["csv", "tsv", "csv.gz", "tsv.gz"])
else:
    vcf_path = st.text_input("VCF/VCF.GZ file path")
    gff_path = st.text_input("GFF3/GFF3.GZ file path")
    traits_path = st.text_input("Traits CSV/TSV file path")

ready = (method == "Upload Files" and vcf_file and gff_file and traits_file) or \
        (method == "Use Local Paths" and vcf_path and gff_path and traits_path)

if ready:
    with st.spinner("🔄 Processing..."):
        try:
            total_start = time.time()

            # Step 1: VCF Parsing
            start = time.time()
            vcf_stream = smart_open(vcf_file if method == "Upload Files" else vcf_path)
            variants_df = filter_variants(vcf_stream, include_intergenic=include_intergenic)
            st.success(f"✅ VCF parsed in {round(time.time() - start, 2)}s - {len(variants_df)} variants")
            st.dataframe(variants_df.head())

            if variants_df.empty:
                st.error("🚫 No variants retained. Try enabling 'Include intergenic variants'.")
                st.stop()

            # Step 2: GFF Parsing
            start = time.time()
            gff_stream = smart_open(gff_file if method == "Upload Files" else gff_path)
            gene_db = build_gene_db(gff_stream)
            st.success(f"✅ GFF parsed in {round(time.time() - start, 2)}s - {len(gene_db)} genes")

            # Step 3: Traits Parsing
            start = time.time()
            traits_stream = smart_open(traits_file if method == "Upload Files" else traits_path)
            traits_df = read_gene_traits(traits_stream)
            st.success(f"✅ Traits parsed in {round(time.time() - start, 2)}s - {len(traits_df)} entries")

            # Step 4: Annotate Genes
            start = time.time()
            annotated_df = annotate_variants_with_genes(variants_df, gene_db, include_intergenic)
            st.success(f"✅ Annotated with genes in {round(time.time() - start, 2)}s")

            # Step 5: Annotate Traits
            start = time.time()
            result_df = annotate_with_traits(annotated_df, traits_df)
            st.success(f"✅ Annotated with traits in {round(time.time() - start, 2)}s")

            # ✅ Final Result
            st.success(f"🎉 Done in {round(time.time() - total_start, 2)} seconds")
            st.dataframe(result_df)

            # 📊 Plot Traits Distribution
            if "Trait" in result_df.columns and not result_df["Trait"].dropna().empty:
                fig, ax = plt.subplots(figsize=(8, 4))
                sns.countplot(data=result_df.dropna(subset=["Trait"]), x="Trait",
                              order=result_df["Trait"].value_counts().index, ax=ax)
                ax.set_title("Variants per Trait")
                ax.set_ylabel("Count")
                ax.set_xlabel("Trait")
                plt.xticks(rotation=45)
                st.pyplot(fig)

            # 💾 Download
            csv = result_df.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Download as CSV", csv, "filtered_variants.csv", "text/csv")

        except Exception as e:
            st.error(f"❌ An error occurred:\n\n{str(e)}")
