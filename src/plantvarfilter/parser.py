# src/plantvarfilter/parser.py

import pandas as pd
import gzip
import io
from typing import Union, TextIO
from pathlib import Path

def smart_open(file_or_path: Union[str, Path, TextIO]) -> TextIO:
    """
    Open a file or stream (optionally gzip-compressed) and return a text stream.
    Supports file-like objects, string paths, and pathlib.Path.
    """
    if isinstance(file_or_path, (str, Path)):
        file_or_path = str(file_or_path)  # Ensure it's a string path
        if file_or_path.endswith(".gz"):
            return gzip.open(file_or_path, "rt", encoding="utf-8")
        else:
            return open(file_or_path, "r", encoding="utf-8")
    elif hasattr(file_or_path, 'read'):
        return file_or_path
    else:
        raise ValueError("Unsupported file input type")



def read_gene_traits(file_or_path: Union[str, TextIO]) -> pd.DataFrame:
    """
    Read a traits file (CSV/TSV) with at least columns: Gene or Gene_ID and Trait.
    """
    stream = smart_open(file_or_path)

    # Try reading as TSV or CSV
    name = getattr(file_or_path, 'name', str(file_or_path)).lower()
    sep = "\t" if name.endswith(".tsv") or name.endswith(".tsv.gz") else ","

    df = pd.read_csv(stream, sep=sep)

    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]

    if "Gene_ID" in df.columns:
        df.rename(columns={"Gene_ID": "Gene"}, inplace=True)

    return df
