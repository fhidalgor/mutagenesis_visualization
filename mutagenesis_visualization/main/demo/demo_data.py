"""
This module contains sample *Screen* objects that can be used to
explore this API.
"""
from typing import Dict
import numpy as np
from pandas import read_pickle
from pandas.core.frame import DataFrame
from mutagenesis_visualization.main.utils.pandas_functions import parse_pivot
from mutagenesis_visualization.main.utils.data_paths import (
    HRAS_RBD_COUNTS_CSV, HRAS_GAPGEF_COUNTS_CSV, DF_BLA_RAW_PKL, DF_SUMO1_RAW_PKL, DF_MAPK1_RAW_PKL,
    DF_UBE2I_RAW_PKL, DF_TAT_PKL, DF_REV_PKL, DF_ASYNUCLEIN_PKL, DF_APH_PKL, DF_B11L5F_RAW_PKL
)


def load_demo_datasets() -> Dict[str, DataFrame]:
    """
    Loads example datasets so the user can play with it.

    Returns
    --------
    data_dict : Dict[str, DataFrame]
        Dictionary that contains the datasets used to create the plots on the documentation.    """
    # Create dictionary where to store data
    data_dict: Dict[str, DataFrame] = {}

    # Retrieve H-Ras datasets and store in dict
    hras_enrichment_rbd = DataFrame(np.genfromtxt(HRAS_RBD_COUNTS_CSV, delimiter=','))
    data_dict['array_hras_RBD'] = hras_enrichment_rbd

    hras_enrichment_gapgef = DataFrame(np.genfromtxt(
        HRAS_GAPGEF_COUNTS_CSV,
        delimiter=',',
    ))
    data_dict['array_hras_gapgef'] = hras_enrichment_gapgef

    # Beta lactamase data
    df_bla_raw = read_pickle(DF_BLA_RAW_PKL)
    data_dict['df_bla'], _ = parse_pivot(df_bla_raw, col_data='DMS_amp_625_(b)')

    # Sumo
    df_sumo1_raw = read_pickle(DF_SUMO1_RAW_PKL)
    data_dict['df_sumo1'], _ = parse_pivot(df_sumo1_raw, col_data='DMS')

    # MAPK1
    df_mapk1_raw = read_pickle(DF_MAPK1_RAW_PKL)
    data_dict['df_mapk1'], _ = parse_pivot(df_mapk1_raw, col_data='DMS_DOX')

    # UBE2I
    df_ube2i_raw = read_pickle(DF_UBE2I_RAW_PKL)
    data_dict['df_ube2i'], _ = parse_pivot(df_ube2i_raw, col_data='DMS')

    # TAT
    data_dict['df_tat'] = read_pickle(DF_TAT_PKL)

    # REV
    data_dict['df_rev'] = read_pickle(DF_REV_PKL)

    # asynuclein
    data_dict['df_asynuclein'] = read_pickle(DF_ASYNUCLEIN_PKL)

    # APH
    data_dict['df_aph'] = read_pickle(DF_APH_PKL)

    # b11L5
    df_b11l5f_raw = read_pickle(DF_B11L5F_RAW_PKL)
    data_dict['df_b11l5f'], _ = parse_pivot(df_b11l5f_raw, col_data='relative_tryp_stability_score')

    return data_dict
