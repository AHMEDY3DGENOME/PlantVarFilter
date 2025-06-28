import gzip
import io
import pandas as pd

def smart_open(file_or_path):
    """
    Open a file-like object or path (compressed or not) and return a text stream.
    Supports Streamlit uploads and local file paths.
    """
    if isinstance(file_or_path, str):
        if file_or_path.endswith(".gz"):
            return io.TextIOWrapper(gzip.open(file_or_path, "rb"), encoding="utf-8")
        else:
            return open(file_or_path, "r", encoding="utf-8")

    elif hasattr(file_or_path, 'read'):
        # For Streamlit uploaded files
        if hasattr(file_or_path, 'name') and file_or_path.name.endswith(".gz"):
            return io.TextIOWrapper(gzip.GzipFile(fileobj=file_or_path), encoding="utf-8")
        else:
            return io.TextIOWrapper(file_or_path, encoding="utf-8")

    else:
        raise ValueError("Unsupported file input type.")

def read_gene_traits(file_or_path):
    """
    Read trait annotation file (CSV, TSV, compressed or not) into a DataFrame.
    The file must contain at least a 'Gene' or 'Gene_ID' column.
    """
    stream = smart_open(file_or_path)

    # Determine extension
    if isinstance(file_or_path, str):
        name = file_or_path
    else:
        name = getattr(file_or_path, 'name', 'unknown.csv')

    # Try to load file with automatic separator detection
    if name.endswith(".tsv") or name.endswith(".tsv.gz"):
        df = pd.read_csv(stream, sep="\t")
    else:
        try:
            df = pd.read_csv(stream, sep=None, engine="python")
        except Exception:
            # fallback to comma
            stream.seek(0)
            df = pd.read_csv(stream, sep=",")

    # Remove duplicate columns if exist
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]

    # Rename Gene_ID to Gene if needed
    if "Gene_ID" in df.columns and "Gene" not in df.columns:
        df.rename(columns={"Gene_ID": "Gene"}, inplace=True)

    # Ensure 'Gene' column exists
    if "Gene" not in df.columns:
        raise ValueError("Trait file must contain 'Gene' or 'Gene_ID' column.")

    return df
