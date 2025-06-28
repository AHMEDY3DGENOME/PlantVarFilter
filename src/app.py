import streamlit as st
import time
import matplotlib.pyplot as plt
import seaborn as sns
import gzip
from filter import improved_filter_variants as filter_variants, build_gene_db, annotate_variants_with_genes, annotate_with_traits
from parser import smart_open, read_gene_traits

st.set_page_config(page_title="PlantVarFilter", layout="wide")
st.title("\U0001F33F PlantVarFilter: Variant Filtering for Plant Genomics")

method = st.radio("Select file input method:", ["Upload Files", "Use Local Paths"])

vcf_file = gff_file = traits_file = None
vcf_path = gff_path = traits_path = ""

include_intergenic = st.checkbox("Include intergenic variants", value=True)

if method == "Upload Files":
    vcf_file = st.file_uploader("Upload VCF or VCF.GZ File", type=["vcf", "vcf.gz"])
    gff_file = st.file_uploader("Upload GFF3 or GFF3.GZ File", type=["gff", "gff3", "gff.gz", "gff3.gz"])
    traits_file = st.file_uploader("Upload Trait CSV/TSV File", type=["csv", "tsv", "csv.gz", "tsv.gz"])
else:
    vcf_path = st.text_input("Enter path to VCF or VCF.GZ file")
    gff_path = st.text_input("Enter path to GFF3 or GFF3.GZ file")
    traits_path = st.text_input("Enter path to traits CSV/TSV file")

if (method == "Upload Files" and vcf_file and gff_file and traits_file) or \
   (method == "Use Local Paths" and vcf_path and gff_path and traits_path):

    with st.spinner("Processing files..."):
        try:
            total_start = time.time()

            # Step 1: Filter VCF variants
            start = time.time()
            if method == "Upload Files":
                vcf_stream = gzip.open(vcf_file) if vcf_file.name.endswith(".gz") else vcf_file
                variants_df = filter_variants(vcf_stream, include_intergenic=include_intergenic)
            else:
                with open(vcf_path, "rb") as f:
                    vcf_stream = gzip.open(f) if vcf_path.endswith(".gz") else f
                    variants_df = filter_variants(vcf_stream, include_intergenic=include_intergenic)
            st.write(f"✅ VCF filtered in {round(time.time() - start, 2)} seconds")
            st.write(f"🔍 Variants after filtering: {len(variants_df)}")
            st.dataframe(variants_df.head())

            # ✅ Stop if no variants found
            if variants_df.empty:
                st.error("🚨 No variants were retained after filtering. Please check your VCF file or try enabling intergenic variants.")
                st.stop()

            # Step 2: Build gene DB from GFF
            start = time.time()
            if method == "Upload Files":
                gff_stream = gzip.open(gff_file) if gff_file.name.endswith(".gz") else gff_file
                gene_db = build_gene_db(gff_stream)
            else:
                with open(gff_path, "rb") as f:
                    gff_stream = gzip.open(f) if gff_path.endswith(".gz") else f
                    gene_db = build_gene_db(gff_stream)
            st.write(f"✅ GFF DB built in {round(time.time() - start, 2)} seconds")

            # Step 3: Load traits
            start = time.time()
            if method == "Upload Files":
                traits_stream = gzip.open(traits_file) if traits_file.name.endswith(".gz") else traits_file
                traits_df = read_gene_traits(traits_stream)
            else:
                if traits_path.endswith(".gz"):
                    with gzip.open(traits_path, "rt") as f:
                        traits_df = read_gene_traits(f)
                else:
                    traits_df = read_gene_traits(traits_path)
            st.write(f"✅ Traits loaded in {round(time.time() - start, 2)} seconds")

            # Step 4: Annotate variants with genes
            start = time.time()
            annotated_df = annotate_variants_with_genes(variants_df, gene_db, include_intergenic=include_intergenic)
            st.write(f"✅ Variants annotated with genes in {round(time.time() - start, 2)} seconds")

            # Step 5: Annotate with traits
            start = time.time()
            result_df = annotate_with_traits(annotated_df, traits_df, keep_unmatched=True)
            st.write(f"✅ Annotation with traits completed in {round(time.time() - start, 2)} seconds")

            # ✅ Final display
            st.success(f"🎉 Done in {round(time.time() - total_start, 2)} seconds")
            st.dataframe(result_df)

            # Optional plot
            if not result_df["Trait"].dropna().empty:
                fig, ax = plt.subplots()
                sns.countplot(data=result_df.dropna(subset=["Trait"]), x="Trait",
                              order=result_df["Trait"].value_counts().index, ax=ax)
                ax.set_title("Variants by Trait")
                ax.set_ylabel("Count")
                plt.xticks(rotation=45)
                st.pyplot(fig)

            # Download CSV
            csv = result_df.to_csv(index=False).encode("utf-8")
            st.download_button("Download Results as CSV", csv, "filtered_variants.csv", "text/csv")

        except Exception as e:
            st.error(f"❌ Error during processing: {e}")
