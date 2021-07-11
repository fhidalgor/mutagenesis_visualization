def plot_miniheatmap(
    self,
    offset=0,
    background_correction=True,
    output_file: Union[None, str, Path] = None,
    **kwargs,
):
    """
    Generate a miniheatmap plot enrichment scores of mutagenesis selection
    assays.

    Parameters
    ----------
    self : object from class *Screen*

    offset : int, default 0
        Will group columns by residues. If the offset is not 0, it will use the values
        of the n+offset to group by. For example, you may want to see what happens when
        you have a Proline in front of the mutated residue. The algorithm can report
        the difference between the calculated value and the mean score for that particular
        substitution.
        Offset of 1 means that you evaluate the effect of following residue n+1 on n.
        Offset of -1 means that you look at the previous residue (n-1 on n).

    background_correction : boolean, default True
        If offset is nonzero, whether subtract the average effect of a substitution or not.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    **kwargs : other keyword arguments
        return_plot_object : boolean, default False
            If true, will return plotting objects (ie. fig, ax).

    Returns
    ----------
    fig, ax, cb : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    """

    # load font parameters
    code_kwargs._font_parameters()

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)

    # load labels
    temp_kwargs['color_sequencelabels'] = _labels(self.start_position)[0]
    temp_kwargs['number_sequencelabels'] = _labels(self.start_position)[1]

    # do offset if appropriate
    dataframe_stopcodons = _transform_dataset_offset(self, offset)

    # calculate condensed heatmap
    dataset = _condense_heatmap(
        dataframe_stopcodons, temp_kwargs['neworder_aminoacids']
    )

    # if offset is not 0
    if background_correction and offset != 0:
        if '*' in temp_kwargs['neworder_aminoacids']:
            temp_kwargs['neworder_aminoacids'].remove('*')
        # do offset, no stop codons
        dataset = _normalize_neighboreffect(
            self, offset, temp_kwargs['neworder_aminoacids']
        )

    fig, ax, cb = _plot_miniheatmap(dataset, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax, cb


def _condense_heatmap(df, new_order):
    '''
    Converts the np.array with stored enrichment scores into the condensed heatmap
    '''
    # Convert dataset to df
    df = df.copy()
    df.drop(['Position'], axis=1, inplace=True)

    # Group by sequence and aminoacid, and then pivot table
    df_grouped = df.groupby(['Sequence', 'Aminoacid'], sort=False).mean()
    df_pivoted = df_grouped.pivot_table(
        values='Score', index='Aminoacid', columns='Sequence'
    )
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
        trimmedsequence = sequence[start_position - 1 + offset:len(dataset[0]) +
                                   start_position - 1 + offset]
    else:
        sequence = 'X' * (np.absolute(offset)) + sequence
        trimmedsequence = sequence[start_position - 1:len(dataset[0]) +
                                   start_position - 1]

    return trimmedsequence


def _transform_dataset_offset(self, offset, stopcodons=True):
    '''
    Generate a dataframe with the sequence offset. Reutilizes _transform_dataset
    '''
    # Add offset sequence
    offset_sequence = _offset_sequence(
        self.dataset, self.sequence_raw, self.start_position, offset
    )
    df = self.dataframe_stopcodons.copy(
    ) if stopcodons is True else self.dataframe.copy()

    # Copy old sequence
    df['Sequence_old'] = df['Sequence']
    # Count amino acids
    aa_number = len(set(df['Aminoacid']))
    # Generate new offset sequence
    df['Sequence'] = np.ravel([[aa] * aa_number for aa in offset_sequence])

    # Drop rows with X
    df.drop(df.index[df['Sequence'] == 'X'], inplace=True)

    return df
