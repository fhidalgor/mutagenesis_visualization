def plot_roc(
    self, df_class=None, mode='pointmutant',output_file: Union[None, str, Path] = None, **kwargs
):
    """
    Generates ROC AUC plot. It compares enrichment scores to some labels that
    the user has specified.

    Parameters
    -----------
    self : object from class *Screen*

    df_class: Pandas dataframe
        A dataframe that contains a column of variants labeled 'Variant' with a column labeled 'Class'
        containing the true class of that mutation. The true class can also be an input when creating the object.

    mode : str, default 'pointmutant'
        Specify what enrichment scores to show. If mode = 'mean', it will show the mean of
        each position. If mode = 'A', it will show the alanine substitution profile. Can be
        used for each amino acid. Use the one-letter code and upper case.

    output_file : str, default None
        If you want to export the generated graph, add the path and name of the file.
        Example: 'path/filename.png' or 'path/filename.svg'.

    Returns
    --------
    fig, ax : matplotlib figure and subplots
        Needs to have return_plot_object==True. By default they do
        not get returned.

    """
    # Chose mode:
    df_grouped = _select_grouping(self.dataframe, mode)

    # Use default class
    if df_class is None:
        df_class = self.roc_df

    # Merge dataframe with classes
    df = _mergeclassvariants(self.dataframe, df_class, mode)

    # Calculate ROC parameters
    fpr, tpr, auc, _ = _rocauc(df)

    # update kwargs
    temp_kwargs = copy.deepcopy(code_kwargs.kwargs())
    temp_kwargs.update(kwargs)
    temp_kwargs['figsize'] = kwargs.get('figsize', (2.5, 2.5))

    # import parameters
    code_kwargs._parameters()

    # create figure
    fig, ax = plt.subplots(figsize=temp_kwargs['figsize'])
    lw = 2
    plt.plot(fpr, tpr, color='k', lw=lw, label='AUC = %0.2f' % auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')

    # Graph limits
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    tick_spacing = 0.2
    ax.xaxis.set_major_locator(
        ticker.MultipleLocator(tick_spacing)
    )  # Plt ticks
    ax.yaxis.set_major_locator(
        ticker.MultipleLocator(tick_spacing)
    )  # Plt ticks

    # Axis labels
    plt.title(temp_kwargs['title'], fontsize=12, fontname='Arial', color='k')
    plt.ylabel(
        'True Positive Rate',
        fontsize=12,
        fontname="Arial",
        color='k',
        labelpad=0
    )
    plt.xlabel('False Positive Rate', fontsize=12, fontname="Arial", color='k')

    # Legend
    plt.legend(
        loc='lower right', handlelength=0, handletextpad=0, frameon=False
    )

    # save file
    code_utils._save_work(fig, output_file, temp_kwargs)

    # return matplotlib object
    if temp_kwargs['return_plot_object']:
        return fig, ax

    if temp_kwargs['show']:
        plt.show()

def _select_grouping(df, mode):
    '''
    Choose the subset of substitutions based on mode input.
    For example, if mode=='A', then return data for Alanine.

    '''
    # convert to upper case
    mode = mode.upper()

    # Select grouping
    if mode =='POINTMUTANT':
        pass
    elif mode == 'MEAN':
        df = df.groupby('Position', as_index=False).mean()
    else:
        df = df.loc[self.dataframe['Aminoacid'] == mode].copy()

    return df

def _rocauc(df):
    '''
    Calculate roc rates and auc.

    The input is a dataframe that contains [Variants,Class,Score]
    '''
    fpr, tpr, thresholds = metrics.roc_curve(
        df['Class'], df['Score'], drop_intermediate=True
    )
    auc = metrics.roc_auc_score(df['Class'], df['Score'])
    return fpr, tpr, auc, thresholds


def _mergeclassvariants(df_score, df_class, mode):
    '''
    Merge the input dataframe containing the class (true score) for variants and the enrichment scores
    '''
    # convert to upper case
    mode = mode.upper()

    if mode == 'POINTMUTANT':
        # Cut other data
        df_class = df_class[['Variant', 'Class']].copy()
        # Merge DMS with true score dataset
        df_merged = pd.merge(df_score, df_class, on=['Variant'], how='left')
    else:
            # Cut other data
        df_class = df_class[['Position','Class']].copy()
        df_merged = pd.merge(df_score, df_class, on=['Position'], how='left')

    # Drop rows with Nan values
    df_merged.dropna(inplace=True)

    return df_merged
