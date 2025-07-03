import pytest
import pandas as pd
import matplotlib.pyplot as plt
from plantvarfilter.visualize import plot_trait_counts

def test_plot_trait_counts_basic():
    data = {
        "Trait": ["A", "B", "A", "C", "B", "B", None]
    }
    df = pd.DataFrame(data)

    fig = plot_trait_counts(df, trait_column="Trait")
    assert fig is not None
    assert isinstance(fig, plt.Figure)

def test_plot_trait_counts_empty_after_dropna():
    df = pd.DataFrame({"Trait": [None, None]})
    with pytest.raises(ValueError, match="No traits to plot."):
        plot_trait_counts(df, trait_column="Trait")

def test_plot_trait_counts_with_custom_ax():
    data = {
        "Trait": ["X", "Y", "X", "Z"]
    }
    df = pd.DataFrame(data)
    fig, ax = plt.subplots()
    returned_fig = plot_trait_counts(df, trait_column="Trait", ax=ax)
    assert returned_fig is fig
