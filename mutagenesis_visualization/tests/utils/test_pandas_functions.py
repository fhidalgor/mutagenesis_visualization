"""
This module tests the utils.
"""
import numpy as np
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.utils.pandas_functions import (
    return_common_elements, process_mean_residue, process_by_pointmutant, color_data
)


def test_return_common_elements() -> None:
    """
    Test for return common elements of two lists function.
    """
    # lists
    list_a = list('abcde')
    list_b = list('efghi')
    list_c = list('jklmn')
    list_d = list('mklj')

    # error message
    error_message = 'return_common_elements not returning the common elements of two lists.'

    # assert
    assert return_common_elements(list_a, list_b) == ['e'], error_message
    assert return_common_elements(list_a, list_c) == [], error_message
    assert return_common_elements(list_c, list_d) == list('jklm'), error_message


def test_process_by_pointmutant() -> None:
    """
    Testing that output type is a dataframe.
    """
    # Create dataframes as attributes of the objects
    dataframe_1: DataFrame = DataFrame(
        np.array([[1, 2], [4, 5], [7, 8]]), columns=['Score_NaN', 'Variant']
    )
    dataframe_2: DataFrame = DataFrame(np.array([[7, 8], [9, 0]]), columns=['Score_NaN', 'Variant'])

    # Call the function we are testing
    df_output: DataFrame = process_by_pointmutant(dataframe_1, dataframe_2)

    # Assert
    assert len(df_output) == 2, 'Truncation of longer dataset is not working properly.'


def test_process_meanresidue() -> None:
    """
    Testing full capabilities of function.
    """
    # Create dataframes as attributes of the objects
    dataframe_1: DataFrame = DataFrame(
        np.array([[1, 2], [1, 6], [2, 8], [2, 4]]), columns=['Position', 'Score']
    )
    dataframe_2: DataFrame = DataFrame(
        np.array([[7, 8], [9, 0], [1, 6]]), columns=['Position', 'Score']
    )
    expected_answer: DataFrame = DataFrame(
        np.array([[4.0, 6.0, int(1), -2.0], [6.0, 8.0, int(2), -2.0]]),
        columns=['dataset_1', 'dataset_2', 'Position', 'd1 - d2']
    ).astype({"Position": int})
    # Call the function we are testing
    df_output: DataFrame = process_mean_residue(dataframe_1, dataframe_2)

    # Assert
    assert df_output.equals(expected_answer), 'error in _process_mean_residue'


def test_color_data() -> None:
    """
    Testing full capabilities of function.
    """
    df_input: DataFrame = DataFrame()
    df_input['Score'] = [1, 2, 3, 0, -1, -2, -3]
    df_input['Expected_Answer'] = ['red'] * 3 + ['blue'] * 4
    df_input['Calculated_answer'] = [
        color_data(df_input.loc[i], 'red', 'blue') for i in range(0, len(df_input['Score']))
    ]
    assert (df_input['Expected_Answer'] == df_input['Calculated_answer']
            ).all(), 'error when assigning a color'
