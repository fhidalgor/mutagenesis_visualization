from mutagenesis_visualization.main.utils.pandas_functions import df_rearrange


def select_grouping(df, mode):
    """
    Choose the subset of substitutions based on mode input.
    For example, if mode=='A', then return data for Alanine.

    """
    # convert to upper case
    mode = mode.upper()

    # Select grouping
    if mode == 'POINTMUTANT':
        pass
    elif mode == 'MEAN':
        df = df.groupby('Position', as_index=False).mean()
    else:
        df = df.loc[self.dataframe['Aminoacid'] == mode].copy()

    return df


def roc_auc(df):
    """
    Calculate roc rates and auc.

    The input is a dataframe that contains [Variants,Class,Score]
    """
    fpr, tpr, thresholds = metrics.roc_curve(df['Class'], df['Score'], drop_intermediate=True)
    auc = metrics.roc_auc_score(df['Class'], df['Score'])
    return fpr, tpr, auc, thresholds


def merge_class_variants(df_score, df_class, mode):
    """
    Merge the input dataframe containing the class (true score) for
    variants and the enrichment scores
    """
    # convert to upper case
    mode = mode.upper()

    if mode == 'POINTMUTANT':
        # Cut other data
        df_class = df_class[['Variant', 'Class']].copy()
        # Merge DMS with true score dataset
        df_merged = pd.merge(df_score, df_class, on=['Variant'], how='left')
    else:
        # Cut other data
        df_class = df_class[['Position', 'Class']].copy()
        df_merged = pd.merge(df_score, df_class, on=['Position'], how='left')

    # Drop rows with Nan values
    df_merged.dropna(inplace=True)

    return df_merged


def condense_heatmap(df, new_order):
    """
    Converts the np.array with stored enrichment scores into the condensed heatmap
    """
    # Convert dataset to df
    df = df.copy()
    df.drop(['Position'], axis=1, inplace=True)

    # Group by sequence and aminoacid, and then pivot table
    df_grouped = df.groupby(['Sequence', 'Aminoacid'], sort=False).mean()
    df_pivoted = df_grouped.pivot_table(values='Score', index='Aminoacid', columns='Sequence')
    df_pivoted.reset_index(drop=False, inplace=True)

    # Sort in y axis desired order
    df_pivoted['Aminoacid'] = pd.Categorical(df_pivoted['Aminoacid'], new_order)
    df_pivoted = df_pivoted.sort_values(by=['Aminoacid'])

    # Sort in x axis desired order
    x_order = code_utils._common(new_order, list(df_pivoted.columns))

    # Drop amino acid column
    data_dropped = df_pivoted.drop(['Aminoacid'], axis=1)

    return data_dropped[x_order]


def _offset_sequence(dataset, sequence, start_position, offset):
    """
    Internal function that offsets the input sequence.

    Parameters
    -----------
    dataset, sequence, start_position, offset

    Returns
    --------
    string containing trimmed sequence

    """
    # Deep copy sequence
    sequence = copy.deepcopy(sequence)

    # truncate sequence
    if offset > 0:
        sequence = sequence + 'X' * np.absolute(offset)
        trimmedsequence = sequence[start_position - 1 + offset : len(dataset[0]) + start_position -
                                   1 + offset]
    else:
        sequence = 'X' * (np.absolute(offset)) + sequence
        trimmedsequence = sequence[start_position - 1 : len(dataset[0]) + start_position - 1]

    return trimmedsequence


def transform_dataset_offset(self, offset, stopcodons=True):
    """
    Generate a dataframe with the sequence offset. Reutilizes _transform_dataset
    """
    # Add offset sequence
    offset_sequence = _offset_sequence(self.dataset, self.sequence_raw, self.start_position, offset)
    df = self.dataframe_stopcodons.copy() if stopcodons is True else self.dataframe.copy()

    # Copy old sequence
    df['Sequence_old'] = df['Sequence']
    # Count amino acids
    aa_number = len(set(df['Aminoacid']))
    # Generate new offset sequence
    df['Sequence'] = np.ravel([[aa] * aa_number for aa in offset_sequence])

    # Drop rows with X
    return df.drop(df.index[df['Sequence'] == 'X'], inplace=True)


def _normalize_neighboreffect(self, offset, neworder):
    """
   For every residue, subtract the average effect of a substitution
   Returns a normalized dataframe
   """
    aalist = list('ACDEFGHIKLMNPQRSTVWY')
    # Add offset sequence to df
    df = _transform_dataset_offset(self, offset, False)

    # calculate mean effect using condensed heatmap
    mean = _condense_heatmap(self.dataframe, aalist)

    df_normalized = pd.DataFrame()
    for aa in aalist:
        # Choose the neighbors of an aa
        aa_neighbors = df.loc[df['Sequence'] == aa]
        # Do the mean substitution of amino acids that are repeated
        aa_neighbors = aa_neighbors.groupby(['Sequence_old', 'Aminoacid'], as_index=False).mean()
        # Make into table
        aa_neighbors_pivoted = aa_neighbors.pivot_table(
            values='Score', index='Aminoacid', columns='Sequence_old'
        )
        aa_neighbors_pivoted.reset_index(drop=True, inplace=True)
        # Get the mean of the amino acids that appear in the aa_neighbors subset
        mean_neighbors = mean[list(aa_neighbors_pivoted.columns)]
        # Subtract average effect and do mean
        df_normalized[aa] = (aa_neighbors_pivoted - mean_neighbors).mean(axis=1)

    # Sort by aa
    df_normalized = df_normalized[neworder]
    # Sort in y axis desired order
    df_normalized = _sort_yaxis_aminoacids(df_normalized, neworder, aalist)
    return df_normalized


def _sort_yaxis_aminoacids(df, neworder, oldorder=list('ACDEFGHIKLMNPQRSTVWY')):
    # Sort in y axis desired order
    df['Aminoacid_new'] = oldorder
    df['Aminoacid_new'] = pd.Categorical(df['Aminoacid_new'], neworder)
    df.sort_values(by=['Aminoacid_new'], inplace=True)
    df.drop(['Aminoacid_new'], inplace=True, axis=1)

    return df
