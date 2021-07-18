"""
This module tests the classes.
"""
import traceback

import pandas as pd
import numpy as np

import mutagenesis_visualization.main.scripts.code_class as code_class

def test_counts()-> None:

    # fake dataframe
    df = pd.DataFrame(
        np.random.rand(21, 10) * 100, index=list('ACDEFGHIKLMNPQRSTVWY*')
    )
    aminoacids = list('ACDEFGHIKLMNPQRSTVWY*')

    # Define aux function
    def _test_Counts_output(parameters)->bool:
        try:
            code_class.Counts(**parameters)
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return True
        return False

#    assert _test_Class_output(df) == False, 'Error when generating class Counts'
    list_params = [{'df': df},
                   {'df': df, 'start_position': 1},
                   {'df': df, 'aminoacids': aminoacids}]
        # Assert
    for parameters in list_params:
        assert _test_Counts_output(
            parameters
        ) == False, "class Counts failed with {} parameters".format(parameters)


def test_screen():

    # fake dataframe
    df = pd.DataFrame(np.random.rand(21, 10))
    sequence = 'MTEYKRVVVLL'
    secondary = ['β1'] * len(sequence)

    # Define aux function
    def _test_Screen_output(parameters):
        try:
            code_class.Screen(**parameters)
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return True
        return False

    list_params = [{'dataset': df, 'sequence': sequence},
                   {'dataset': df, 'sequence': sequence}]

    # Assert
    for parameters in list_params:
        assert _test_Screen_output(
            parameters
        ) == False, "class Screen failed with {} parameters".format(parameters)
