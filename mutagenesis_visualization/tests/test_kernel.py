"""
This module tests the kernel and histogram methods.
"""
import traceback
import logging
from typing import List
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects

DEMO_OBJECTS: DemoObjects = DemoObjects()
OBJ_TEST = DEMO_OBJECTS.hras_rbd


def test_kernel() -> None:
    """
    Test the kernel method.
    """
    def _test_kernel_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.kernel(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [
        {'show': False},
        {'cumulative': True, 'show': False, 'close': True},
        {
            'y_label': 'testing', 'title': 'this is a test', 'x_label': 'testing', 'xscale': 2,
            'show': False, 'close': True
        }  #y_label does not change
    ]

    # Assert
    for parameters in list_params:
        assert _test_kernel_output(
            parameters
        ) is False, "plot_kernel failed with {} parameters".format(parameters)


def test_histogram() -> None:
    """
    Test the histogram method.
    """
    def _test_histogram_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.histogram(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [
        {'show': False},
        {'population': 'SNV', 'show': False, 'close': True},
        {'population': 'nonSNV', 'show': False, 'close': True},
    ]

    # Assert
    for parameters in list_params:
        assert _test_histogram_output(
            parameters
        ) == False, "plot_histogram failed with {} parameters".format(parameters)


def test_multiple_kernel() -> None:
    """
    Test of the multiple kernels method.
    """
    def _test_multiple_kernel(parameters: dict) -> bool:
        try:
            OBJ_TEST.multiple_kernel(DEMO_OBJECTS.bla, **parameters)
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [
        {
            'label_kernels': ["H-Ras", "BLA"],
            'show': False,
            'close': True,
        },
        {
            'label_kernels': ["H-Ras", "BLA"], 'figsize': (3, 2.5), 'y_label': r'$∆E^i_x$', 'show':
            False, 'close': True, 'title': 'go bears'
        },
    ]

    # Assert
    for parameters in list_params:  # Loop over the parameters
        assert _test_multiple_kernel( # Assert that that set of parameters works on that object
            parameters,
        ) is False, "plot_multiplekernel failed with {} parameters".format(
                parameters
            )
