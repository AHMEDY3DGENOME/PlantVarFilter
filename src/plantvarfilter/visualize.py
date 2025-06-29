import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from typing import Optional

def plot_trait_counts(df: pd.DataFrame, trait_column: str = "Trait", ax: Optional[plt.Axes] = None) -> plt.Figure:
    """
    Plot the count of variants per trait.

    Parameters:
    - df: DataFrame that includes trait annotations.
    - trait_column: column name in the dataframe that holds the trait info.
    - ax: optional matplotlib Axes to plot into.

    Returns:
    - fig: the matplotlib Figure object.
    """
    df_clean = df.dropna(subset=[trait_column])
    if df_clean.empty:
        raise ValueError("No traits to plot.")

    if ax is None:
        fig, ax = plt.subplots(figsize=(10, 5))
    else:
        fig = ax.figure

    sns.countplot(
        data=df_clean,
        x=trait_column,
        order=df_clean[trait_column].value_counts().index,
        ax=ax
    )

    ax.set_title("Variants by Trait")
    ax.set_ylabel("Count")
    ax.set_xlabel("Trait")
    plt.xticks(rotation=45)
    fig.tight_layout()
    return fig
