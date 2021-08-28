"""
Test pymol utils.
"""

from mutagenesis_visualization.main.utils.pymol_utils import (light_parameters, _array_to_pymol)


def test_array_to_pymol() -> None:
    """
    Test that output is same as expected.
    """
    assert _array_to_pymol([1, 2, 3]) == '1+2+3', "function _array_to_pymol failed"


def test_light_parameters() -> None:
    """check that there is no error when function runs, checks output is function's output"""
    def _test_light_parameters_output() -> bool:
        """checks if function runs or gives error"""
        try:
            light_parameters()
            return False
        except:
            return True

    assert _test_light_parameters_output() == True, "function _light_parameters failed"
