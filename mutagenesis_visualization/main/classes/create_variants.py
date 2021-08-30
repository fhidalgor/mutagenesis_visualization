"""
This module has the class to generate the variants to order to Twist
or other DNA synthesis companies.
"""
from pathlib import Path
from typing import Union, List
from pandas.core.frame import DataFrame


class CreateVariants:
    """
    Class to create variants for DNA synthesis.
    """
    def __init__(self) -> None:
        """
        Start.
        """
        self.dna: str = ""
        self.codon_list: List[str] = []
        self.seq_list: List[str] = []
        self.df_output: DataFrame = DataFrame()

    def __call__(self, dna: str, codon_list: Union[list, str]) -> DataFrame:
        """
        Generate a list of all point mutants given a dna sequence and a list
        of codons.

        Parameters
        -----------
        dna : str,
            Contains the DNA sequence of the allele of reference (usually wild-type).

        codon_list : list or str
            Input a list of the codons that were used to create point
            mutations. Example: ["GCC", "GCG", "TGC"]. It is important to
            know that the order of the codon_list will determine the output order.

        Returns
        --------
        df_output : pandas dataframe
            Dataframe containing the generated sequences.

        """
        # Make upper case in case input was lower case
        self.dna = dna.upper()
        self.codon_list = [item.upper() for item in codon_list]

        # Generate list of variants
        self.seq_list = self._enumerate_variants_2()

        # Make dataframe
        self.df_output = DataFrame()
        self.df_output['Sequences'] = self.seq_list
        return self.df_output

    def export_file(self, output_file: Union[str, Path]) -> None:
        """
        Parameters
        -----------
        output_file : str
            File where the list of primers will be exported to. Only exports
            to 'xlsx', 'fasta', 'txt'. Example: 'path/sequences.xlsx'.
        """
        if Path(output_file).suffix == '.xlsx':
            self.df_output.to_excel(Path(output_file), sheet_name='Variants', index=False)
        elif Path(output_file).suffix == '.fasta' or Path(output_file).suffix == '.txt':
            self._list_to_fasta(output_file)

    def _enumerate_variants_2(self) -> List[str]:
        """
        Copy of _enumerate variants with slight changes. Does return the
        wild-type sequence as the first item of the list.
        """
        # Create list with codons of sequence
        wtseq: List[str] = [self.dna[i : i + 3] for i in range(0, len(self.dna), 3)]

        # List of sequences
        seq_list: List[str] = [self.dna]

        # Loop over the dna sequence
        for position in range(0, len(wtseq)):
            for codon in self.codon_list:
                variant = ''.join(wtseq[0 : position]) + codon + ''.join(wtseq[position + 1 :])
                if variant != self.dna:
                    seq_list.append(variant)
        return seq_list

    def _list_to_fasta(self, output_file: Union[str, Path]) -> None:
        """
        Export list to fasta format.
        """
        # Open file
        with open(str(output_file), "w", encoding="utf-8") as ofile:
            # Loop through list and write into file
            for i, seq in enumerate(self.seq_list):
                line = (">{}\n{}\n").format(str(i), seq)
                ofile.write(line)
