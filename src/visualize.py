import seaborn as sns
import matplotlib.pyplot as plt

def plot_variant_distribution(df, by='Trait'):
    """Plot histogram of variant counts by trait or chromosome."""
    plt.figure(figsize=(10, 6))
    sns.countplot(data=df, x=by, order=df[by].value_counts().index)
    plt.xticks(rotation=45)
    plt.title(f'Variant distribution by {by}')
    plt.tight_layout()
    plt.savefig('variant_distribution.png')
    plt.close()
