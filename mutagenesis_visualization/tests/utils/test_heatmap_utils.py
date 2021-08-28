"""
This module will test the heatmap utils.
"""
from pandas.core.frame import DataFrame
from mutagenesis_visualization.main.utils.heatmap_utils import (
    hierarchical_sort, _helix, labels, _sheet, _loop
)


def test_hierarchical_sort() -> None:
    """
    Test the hierarchical sorting function.
    """
    df_input: DataFrame = DataFrame([[1, 7, 6, 2], [0, 0, 0, 0], [10, 10, 10, 10], [1, 1, 1, 1]])
    result = hierarchical_sort(df_input.T)
    assert (result == [2, 0, 1, 3]).all(), 'columns are not properly sorted out'


def test_helix() -> None:
    """
    Testing function produces matplotlib object.
    """
    assert (
        str(type(_helix(0, 5))) == "<class 'matplotlib.patches.Rectangle'>"
    ), "function _helix failed"


def test_labels() -> None:
    """testing function produces tuple"""
    assert (str(type(labels(1)))) == "<class 'tuple'>", "function _labels failed"


def test_sheet() -> None:
    """
    Testing function prouduces matplotlib object.
    """
    assert (
        str(type(_sheet(1, 5)))
    ) == "<class 'matplotlib.patches.FancyArrow'>", "function _sheet failed"


def test_loop() -> None:
    """
    Testing function produces matplotlib object.
    """
    assert (
        str(type(_loop(1, 5)))
    ) == "<class 'matplotlib.patches.Rectangle'>", "function _loop failed"
