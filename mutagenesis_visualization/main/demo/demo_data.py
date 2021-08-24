"""
This module contains sample *Screen* objects that can be used to
explore this API.
"""
from typing import Dict, Union
import numpy as np
import pandas as pd
from pandas.core.frame import DataFrame
from mutagenesis_visualization.main.utils.pandas_functions import parse_pivot

PDB_5P21: str = "mutagenesis_visualization/data/5p21.pdb"
PDB_1ERM: str = "mutagenesis_visualization/data/1erm.pdb"
PDB_1A5R: str = "mutagenesis_visualization/data/1a5r.pdb"
PDB_1ND4: str = "mutagenesis_visualization/data/1nd4.pdb"
DEMO_FASTA: str = "mutagenesis_visualization/data/Ras_family_trimmed.fasta"
HRAS_FASTQ: str = "mutagenesis_visualization/data/hras.trimmed.fastq"
HRAS_RBD_COUNTS: str = "mutagenesis_visualization/data/hrasRBD_counts.xlsx"
HRAS_GAPGEF_COUNTS: str = "mutagenesis_visualization/data/hrasGAPGEF_counts.xlsx"


def load_demo_datasets() -> Dict[str, DataFrame]:
    """
    Loads example datasets so the user can play with it.

    Returns
    --------
    data_dict : Dict[str, DataFrame]
        Dictionary that contains the datasets used to create the plots on the documentation.

    """
    # Create dictionary where to store data
    data_dict: Dict[str, DataFrame] = {}

    # Retrieve H-Ras datasets and store in dict
    hras_enrichment_rbd = pd.DataFrame(
        np.genfromtxt(
            "mutagenesis_visualization/data/HRas166_RBD.csv",
            delimiter=',',
        )
    )
    data_dict['array_hras_RBD'] = hras_enrichment_rbd

    hras_enrichment_gapgef = pd.DataFrame(
        np.genfromtxt(
            "mutagenesis_visualization/data/HRas166_GAPGEF.csv",
            delimiter=',',
        )
    )
    data_dict['array_hras_gapgef'] = hras_enrichment_gapgef

    # Beta lactamase data
    df_bla_raw = pd.read_pickle("mutagenesis_visualization/data/df_bla_raw.pkl")
    data_dict['df_bla'], _ = parse_pivot(df_bla_raw, col_data='DMS_amp_625_(b)')

    # Sumo
    df_sumo1_raw = pd.read_pickle("mutagenesis_visualization/data/df_sumo1_raw.pkl")
    data_dict['df_sumo1'], _ = parse_pivot(df_sumo1_raw, col_data='DMS')

    # MAPK1
    df_mapk1_raw = pd.read_pickle("mutagenesis_visualization/data/df_mapk1_raw.pkl")
    data_dict['df_mapk1'], _ = parse_pivot(df_mapk1_raw, col_data='DMS_DOX')

    # UBE2I
    df_ube2i_raw = pd.read_pickle("mutagenesis_visualization/data/df_ube2i_raw.pkl")
    data_dict['df_ube2i'], _ = parse_pivot(df_ube2i_raw, col_data='DMS')

    # TAT
    data_dict['df_tat'] = pd.read_pickle("mutagenesis_visualization/data/df_tat.pkl")

    # REV
    data_dict['df_rev'] = pd.read_pickle("mutagenesis_visualization/data/df_rev.pkl")

    # asynuclein
    data_dict['df_asynuclein'] = pd.read_pickle("mutagenesis_visualization/data/df_asynuclein.pkl")

    # APH
    data_dict['df_aph'] = pd.read_pickle("mutagenesis_visualization/data/df_aph.pkl")

    # b11L5
    df_b11l5f_raw = pd.read_pickle("mutagenesis_visualization/data/df_b11l5f_raw.pkl")
    data_dict['df_b11l5f'], _ = parse_pivot(df_b11l5f_raw, col_data='relative_tryp_stability_score')

    return data_dict
