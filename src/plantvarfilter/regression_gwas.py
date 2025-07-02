
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline
import os

def run_regression_gwas(variant_df: pd.DataFrame, traits_df: pd.DataFrame, output_path: str):
    if "Gene" not in variant_df.columns or "Trait_Score" not in traits_df.columns:
        raise ValueError("Required columns 'Gene' or 'Trait_Score' are missing.")

    gene_flags = variant_df["Gene"].value_counts().to_frame("Has_Variant").reset_index()
    gene_flags.rename(columns={"index": "Gene"}, inplace=True)
    gene_flags["Has_Variant"] = 1

    data = traits_df.merge(gene_flags, on="Gene", how="left").fillna({"Has_Variant": 0})
    data["Has_Variant"] = data["Has_Variant"].astype(int)

    feature_cols = ["Has_Variant", "Age", "Environment"]
    X = data[feature_cols]
    y = data["Trait_Score"]

    preprocessor = ColumnTransformer(
        transformers=[("cat", OneHotEncoder(drop="first"), ["Environment"])],
        remainder="passthrough"
    )

    pipeline = make_pipeline(preprocessor, LinearRegression())
    pipeline.fit(X, y)
    coefs = pipeline.named_steps["linearregression"].coef_

    encoded_columns = pipeline.named_steps["columntransformer"].get_feature_names_out()
    results = pd.DataFrame({
        "Feature": encoded_columns,
        "Coefficient": coefs
    })

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    results.to_csv(output_path, index=False)
    print(f"✅ Regression GWAS saved to {output_path}")
