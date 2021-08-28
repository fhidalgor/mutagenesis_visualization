"""
Module that contains the tests for plotly utils.
"""
from pandas.core.frame import DataFrame
from mutagenesis_visualization.main.utils.plotly_utils import (centroid, color_3d_scatter)


def test_centroid() -> None:
    """
    Tests the mean calculated in each of the three dimensions is correct.
    """
    df_test: DataFrame = DataFrame({"x": [1, 2, 9], "y": [4, 5, 9], "z": [5, 8, 8]})
    assert centroid(df_test) == (4.0, 6.0, 7.0), "function _centroid failed"


def test_color_3d_scatter() -> None:
    """
    Check that output is expected: returns new dataframe with Color
    column added.
    """
    df_input: DataFrame = DataFrame({
        "Position": [1, 1, 2, 2, 9, 9], "Aminoacid": ["A", "B", "A", "B", "A", "B"], "Score":
        [5.0, 5.0, -8.0, -8.0, 8.0, 8.0]
    })
    df_calculated: DataFrame = color_3d_scatter(df_input, 'mean', 1, -1)

    df_solution: DataFrame = DataFrame({
        "Position": [1, 2, 9], "Score": [5.0, -8.0, 8.0], "Color": ["red", "blue", "red"]
    })

    assert df_solution.equals(df_calculated) is True

    df_calculated_2 = color_3d_scatter(df_input, 'A', 1, -1)

    df_solution_2 = DataFrame({
        "Position": [1, 2, 9], "Aminoacid": ["A", "A", "A"], "Score": [5.0, -8.0, 8.0], "Color":
        ["red", "blue", "red"]
    })

    assert list(df_solution_2['Color']) == list(df_calculated_2['Color'])
