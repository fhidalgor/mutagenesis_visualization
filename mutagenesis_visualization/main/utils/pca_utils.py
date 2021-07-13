import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA


def auto_text(x, y, textlabels):
    """
    Auto anotates text labels.
    """
    texts = [
        plt.annotate(
            textlabels[i],  # this is the text
            (x[i], y[i]),  # this is the point to label
            textcoords="offset points",  # how to position the text
            xytext=(2, 2),  # distance from text to points (x,y)
            fontsize=8,
            ha='center'
        )  # horizontal alignment can be left, right or center
        for i in range(len(textlabels))
    ]
    return texts


def calculate_clusters(dataset, dimensions, random_state):
    """
    Input the dataframe that needs to be correlated, the dimensions,
    and will calculate PCA descomposition.
    """

    # call pca model
    pca = PCA(n_components=6, random_state=random_state)

    # fit model to df. use aux function correlation_aminoacids
    model = pca.fit(dataset)

    # create df with PCA data
    df_aa = pd.DataFrame((model.components_).T, columns=['PCA1', 'PCA2', 'PCA3', 'PCA4', 'PCA5', 'PCA6'])

    # use kmeans to cluster the two dimensions and color
    dimensionstoplot = df_aa.iloc[:, np.r_[dimensions[0], dimensions[1]]]

    return dimensionstoplot, pca.explained_variance_ratio_


def _grou_by_secondary(df, secondary):
    """
    Groups each secondary motif and makes the mean.

    Returns dataframe. Returns copy
    """
    df = df.copy()
    df.insert(4, 'Secondary', secondary)
    df = df.groupby(['Secondary', 'Aminoacid'], as_index=False).mean()
    df = df.loc[df['Secondary'].str.startswith(('β', 'α'))]
    return df


def calculate_correlation_by_secondary(df, secondary):
    dataset = _grou_by_secondary(df, secondary)
    dataset = dataset.pivot_table(values='Score', index='Secondary', columns='Aminoacid')
    dataset = dataset.T.corr()

    return dataset


def calculate_correlation(df, order_aminoacids):

    dataset = df.copy()
    dataset = dataset.pivot_table(values='Score', index='Position', columns='Aminoacid')
    dataset = dataset.corr()
    dataset = dataset.reindex(index=order_aminoacids)[order_aminoacids]

    return dataset


def calculate_correlation_by_residue(df):

    dataset = df.copy()
    dataset = dataset.pivot_table(values='Score', index='Position', columns='Aminoacid')
    dataset = dataset.T.corr()

    return dataset

def calculate_substitution_correlations(self, aminoacids, groups):
    """
    If a set of residues was chosen, how well would they represent
    the entire population.
    """

    # Get correlation values
    corr_values = _calculate_correlation(self.dataframe, aminoacids)**2
    corr_values.reset_index(inplace=True)

    # Get combinations
    replacement_combinations = list(itertools.product(*groups))

    # Retrieve Correlation values
    df = pd.DataFrame()
    df['Aminoacids'] = list(itertools.chain.from_iterable(groups))
    for combination in replacement_combinations:  # Iterate over a combination
        temp_list = []

        # Iterate over a group of the combination
        for group, aa_selected in zip(groups, combination):
            for aa_nonselected in group:  # Find correlation values from correlation plot
                if aa_nonselected == aa_selected:
                    temp_list.append(1)
                else:
                    temp_list.append(_find_correlation(aa_selected, aa_nonselected, corr_values))
        df[combination] = temp_list  # Store in df
    return _polishdf(df)


def _polishdf(df):
    df_mean = df.copy()
    df_mean = df.mean().to_frame()
    df_mean.reset_index(drop=False, inplace=True)
    df_mean.rename(columns={0: 'R2'}, inplace=True)
    df_mean['Combinations'] = list(df_mean['index'].apply(lambda x: ''.join(x)))
    df_mean.drop(columns=['index'], inplace=True)
    return df_mean


def _find_correlation(aa1, aa2, corr_values):
    return float(corr_values[aa1].loc[corr_values['Aminoacid'] == aa2])
