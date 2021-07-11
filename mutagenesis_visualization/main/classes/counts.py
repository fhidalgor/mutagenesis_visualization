"""

"""

class Counts:
    '''
    *Counts* represents the output of reading a fastq file.

    Parameters
    ----------
    df : pandas dataframe
        Dataframe containing the counts per codon.

    start_position : int, default None
        First position in the protein sequence that will be used for the first column of the
        array. If a protein has been mutated only from residue 100-150, then if start_position = 100,
        the algorithm will trim the first 99 amino acids in the input sequence. The last
        residue will be calculated based on the length of the input array. We have set the default value to 2
        because normally the Methionine in position 1 is not mutated.

    aminoacids : list, default None
        List of aminoacids (in order). Stop codon needs to be '*'.
        If none, it will use the index of the dataframe

    '''
    def __init__(self, df, start_position=None, aminoacids=None):
        self.dataframe = df

        if start_position:
            self.start_position = start_position
            self.positions = np.arange(start_position, len(df.columns) + 1)
        else:  # if none, use the columns of the dataframe
            self.positions = list(df.columns)
            self.start_position = self.positions[0]

        if aminoacids:
            self.aminoacids = aminoacids
        else:  # if aminoacids is none, use the index of the dataframe
            if _is_DNA(df):
                self.aminoacids = _translate_codons(df)
            else:
                self.aminoacids = list(df.index)
        # bar plots

    mean_counts = plot_meancounts
    library_representation = plot_library_representation
