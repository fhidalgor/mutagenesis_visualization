"""
This module has the class to generate primers for saturation mutagenesis.
"""
from typing import Tuple, Union, Optional, List
from pathlib import Path
from pandas.core.frame import DataFrame
from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt


def _reverse_complement(dna: str) -> str:
    """
    Aux function that uses biopython to calculate the reverse
    complement of a DNA string. Includes mixed-base code. More info in
    https://biopython.org/docs/1.75/api/Bio.Seq.html
    """

    # Needs to be converted to str
    return str(Seq(dna).reverse_complement())


def _primer_design(
    dna: str,
    codon: str,
    codon_position: int,
    length_primer: int,
    melting_temp: float,
) -> Tuple[str, str]:
    """
    Aux function to design the degenerate primers given a sequence
    and a codon position. The length of the primer is fixed.

    Returns
    ---------
    forward_primer, reverse_primer
    """
    if not melting_temp:
        forward_primer: str = dna[(codon_position - length_primer): codon_position] + codon + dna[
            (codon_position + 3):(codon_position + length_primer + 3)]
        return forward_primer, _reverse_complement(forward_primer)

    # loop until melting_temp is achieved
    step_size: int = 6
    forward_primer = "AAA"
    while mt.Tm_NN(forward_primer) < melting_temp:
        forward_primer = dna[(codon_position - step_size): codon_position] + codon + dna[
            (codon_position + 3):(codon_position + step_size)]
        step_size += 1
    return forward_primer, _reverse_complement(forward_primer)


def _create_primers_list(
    dna: str,
    start_codon: int,
    end_codon: int,
    codon: str,
    length_primer: int,
    melting_temp: float,
) -> Tuple[List[str], List[str]]:
    """
    Aux function to create list with fp and list with rp.
    """
    forward_primers: List[str] = []
    reverse_primers: List[str] = []
    for codon_position in range(start_codon, end_codon, 3):
        # Create fp, rp for that position
        forward_primer, reverse_primer = _primer_design(
            dna,
            codon,
            codon_position,
            length_primer,
            melting_temp,
        )
        # Append to list
        forward_primers.append(forward_primer)
        reverse_primers.append(reverse_primer)
    return forward_primers, reverse_primers


class GeneratePrimers:
    """
    Class that will generate primers for saturation mutagenesis.
    """
    def __init__(self, dna: str, start: str, end: str) -> None:
        """
        Parameters
        -----------
        dna : string
            DNA sequence containing the protein of study.
            The DNA sequence should also contain at least 15 base pairs
            before the starting ATG and 15 base pairs after the stop codon.

        start : string
            Beginning of the DNA sequence that will be mutageneized.
            For example, if you will start mutating the first methionine,
            copy a few more dna bases so the algorithm can identify it 'ATGACCAGC'.

        end : string
            The algorithm will stop creating primers once it reaches that base.
            For example, if you will stop mutating at the stop codon, copy a
            few more dna bases ie. 'TAAATGATT'.
        """
        self.dna: str = dna.upper()
        self.start: str = start.upper()
        self.end: str = end.upper()
        self.start_codon: int = self.dna.find(self.start)
        self.end_codon: int = self.dna.find(self.end)
        self.df_primers: DataFrame = DataFrame()

    def __call__(
        self,
        codon: str = 'NNS',
        length_primer: int = 15,
        melting_temp: Optional[float] = None,
    ) -> DataFrame:
        """
        Generate primers for saturation mutagenesis.

        Parameters
        -----------
        codon : str, default 'NNS'
            Degenerate codon that will be used to create the primers. Check
            idt's website for a list of all mixed bases and letter code
            (https://www.idtdna.com/pages/products/custom-dna-rna/mixed-bases).
            This parameter should contain 3 letters, although can contain more.

        length_primer: int, default 15
            Number of bases that the primers will have to each side of the
            mutated codon.  Total primer length will be 2*length_primer+3.

        melting_temp : int, default None
            Melting temperature in Celsius of the primers. Will override
            length_primer. If none, primers will have a total length of
            2*length_primer+3

        Returns
        --------
        df : pandas dataframe
            Dataframe containing the primers.

        """

        # Loop through DNA and make a list with fp and second list with rp
        label_fp = ['fp ' + str(i) for i in range(0, int((self.end_codon - self.start_codon) / 3))]
        label_rp = ['rp ' + str(i) for i in range(0, int((self.end_codon - self.start_codon) / 3))]
        forward_primers, reverse_primers = _create_primers_list(
            self.dna, self.start_codon, self.end_codon, codon, length_primer, melting_temp
        )

        # Create dataframe
        self.df_primers = DataFrame({
            'FP_label': label_fp,
            'FP_seq': forward_primers,
            'RP_label': label_rp,
            'RP_seq': reverse_primers,
        })

        return self.df_primers

    def export_file(self, output_file: Union[str, Path]) -> None:
        """
        Parameters
        -----------
        output_file : str
            File where the list of primers will be exported to. Only exports
            to excel or csv. Example: 'path/primers.xlsx'.
        """
        assert Path(output_file).suffix == "csv" or Path(output_file).suffix == "xlsx"
        if Path(output_file).suffix == "csv":
            self.df_primers.to_csv(Path(output_file), index=False)
        else:
            self.df_primers.to_excel(Path(output_file), sheet_name='Primers', index=False)
