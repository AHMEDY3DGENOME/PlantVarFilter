import logging
import os
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import make_pipeline

def _guess_gene_col(df: pd.DataFrame) -> str:
    candidates = [c for c in df.columns if c.lower() in {"gene", "gene_id", "gene_name"}]
    if candidates:
        return candidates[0]
    obj_cols = [c for c in df.columns if pd.api.types.is_object_dtype(df[c])]
    if obj_cols:
        return obj_cols[0]
    raise ValueError("Could not infer gene column. Please include a string column for genes (e.g., 'Gene').")

def run_regression_gwas_dynamic(
    variant_df: pd.DataFrame,
    traits_df: pd.DataFrame,
    output_path: str,
    max_cat_levels: int = 20
):
    var_gene_col = _guess_gene_col(variant_df)
    tr_gene_col  = _guess_gene_col(traits_df)

    variant_genes = set(variant_df[var_gene_col].dropna().astype(str).unique())
    data = traits_df.copy()
    data[tr_gene_col] = data[tr_gene_col].astype(str)
    data["Has_Variant"] = data[tr_gene_col].isin(variant_genes).astype(int)

    numeric_cols_all = [c for c in data.columns if pd.api.types.is_numeric_dtype(data[c])]
    target_traits = [c for c in numeric_cols_all if c != "Has_Variant"]
    if not target_traits:
        logging.warning("No numeric traits found.")
        return

    candidate_features = [c for c in data.columns if c != tr_gene_col]

    def split_features(df, exclude):
        feats = [c for c in candidate_features if c not in exclude]
        num_feats = [c for c in feats if pd.api.types.is_numeric_dtype(df[c])]
        cat_feats = [c for c in feats if c not in num_feats]
        small_cat = []
        for c in cat_feats:
            nuniq = df[c].astype(str).nunique(dropna=True)
            if nuniq <= max_cat_levels:
                small_cat.append(c)
            else:
                logging.debug(f"Skipping high-cardinality categorical column '{c}' ({nuniq} levels)")
        return num_feats, small_cat

    results_rows = []
    coef_rows = []

    for target in target_traits:
        exclude_now = {target, tr_gene_col}
        num_feats, small_cat = split_features(data, exclude_now)

        if "Has_Variant" not in num_feats:
            num_feats = ["Has_Variant"] + num_feats
        else:
            num_feats = ["Has_Variant"] + [c for c in num_feats if c != "Has_Variant"]

        if small_cat:
            preprocessor = ColumnTransformer(
                transformers=[
                    ("cat", OneHotEncoder(drop="first", handle_unknown="ignore"), small_cat)
                ],
                remainder="passthrough"
            )
            X_cols_order = small_cat + num_feats
        else:
            preprocessor = "passthrough"
            X_cols_order = num_feats

        X = data[X_cols_order].copy()
        for c in X.columns:
            if c in num_feats:
                X[c] = pd.to_numeric(X[c], errors="coerce")
        X = X.fillna(0)

        y = pd.to_numeric(data[target], errors="coerce").fillna(0)

        pipe = make_pipeline(preprocessor, LinearRegression())
        try:
            pipe.fit(X, y)
            r2 = pipe.score(X, y)

            if small_cat:
                feature_names = pipe.named_steps["columntransformer"].get_feature_names_out()
            else:
                feature_names = np.array(X_cols_order)

            coefs = pipe.named_steps["linearregression"].coef_
            coef_map = dict(zip(feature_names, coefs))

            has_variant_coef = None
            for fname in feature_names:
                if fname.split("__")[-1] == "Has_Variant" or fname == "Has_Variant":
                    has_variant_coef = coef_map[fname]
                    break

            results_rows.append({
                "Trait": target,
                "R_squared": r2,
                "Coef_Has_Variant": has_variant_coef
            })

            row = {"Trait": target, **{f"Coef[{k}]": v for k, v in coef_map.items()}}
            coef_rows.append(row)

        except Exception as e:
            logging.warning(f"Regression failed for trait '{target}': {e}")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    summary_path = output_path
    details_path = os.path.join(os.path.dirname(output_path), "gwas_regression_coefficients.csv")

    pd.DataFrame(results_rows).sort_values("R_squared", ascending=False).to_csv(summary_path, index=False)
    pd.DataFrame(coef_rows).to_csv(details_path, index=False)

    logging.info(f"Dynamic regression GWAS saved:\n- Summary: {summary_path}\n- Coeffs:  {details_path}")

def run_regression_gwas(variant_df, traits_df, output_path, **kwargs):
    return run_regression_gwas_dynamic(variant_df, traits_df, output_path, **kwargs)
