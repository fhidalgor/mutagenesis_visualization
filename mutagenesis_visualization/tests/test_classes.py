"""
This module tests the classes.
"""
from typing import List
import traceback
import logging
import numpy as np
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.classes.counts import Counts
from mutagenesis_visualization.main.classes.screen import Screen


def test_counts() -> None:
    """
    Test the class *Counts*.
    """
    df_input: DataFrame = DataFrame(
        np.random.rand(21, 10) * 100, index=list('ACDEFGHIKLMNPQRSTVWY*')
    )
    aminoacids: List[str] = list('ACDEFGHIKLMNPQRSTVWY*')

    def _test_counts_output(parameters: dict) -> bool:
        try:
            Counts(**parameters)
        except Exception as e:  # pylint: disable=broad-except
            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    list_params: List[dict] = [{'dataframes':
                                df_input}, {'dataframes': df_input, 'start_position': 1},
                               {'dataframes': df_input, 'aminoacids': aminoacids}]
    for parameters in list_params:
        assert _test_counts_output(
            parameters
        ) is False, "class Counts failed with {} parameters".format(parameters)


def test_screen() -> None:
    """
    Test the class *Screen*.
    """
    df_input: DataFrame = DataFrame(np.random.rand(21, 10))
    sequence: str = 'MTEYKRVVVLL'
    secondary: list = ['β1'] * len(sequence)
    aminoacids: List[str] = list('ACDEFGHIKLMNPQRSTVWY*')

    def _test_screen_output(parameters: dict) -> bool:
        try:
            Screen(**parameters)
        except Exception as e:  # pylint: disable=broad-except
            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    list_params: List[dict] = [{
        'datasets': df_input, 'sequence': sequence, 'aminoacids': aminoacids
    }, {
        'datasets': df_input, 'sequence': sequence, 'aminoacids': aminoacids, 'secondary': secondary
    }]

    for parameters in list_params:
        assert _test_screen_output(
            parameters
        ) is False, "class Screen failed with {} parameters".format(parameters)
