"""
This module will test the heatmap plotting methods.
"""

import traceback
import logging
from typing import List
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects

DEMO_OBJECTS: DemoObjects = DemoObjects()
OBJ_TEST = DEMO_OBJECTS.hras_rbd


def test_heatmap() -> None:
    """
    Test the heatmap plot method.
    """
    def _test_heatmap_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.heatmap(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except
            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [
        {'show': False, 'close': True},
        {'mask_selfsubstitutions': True, 'hierarchical': True, 'show': False, 'close': True},
        {'show_snv': True, 'show': False, 'close': True},
        {'show_cartoon': True, 'show': False, 'close': True},
    ]

    # Assert
    for parameters in list_params:
        assert _test_heatmap_output(
            parameters
        ) is False, "plot_heatmap failed with {} parameters".format(parameters)


def test_heatmap_rows() -> None:
    """
    Test the heatmap rows plot method.
    """
    def _test_heatmap_rows_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.heatmap_rows(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [{"selection": "mean", 'show': False, 'close': True}, {
        "selection": None, 'mask_selfsubstitutions': True, 'show': False, 'close': True
    }]

    # Assert
    for parameters in list_params:
        assert _test_heatmap_rows_output(
            parameters
        ) is False, "plot_heatmap_rows failed with {} parameters".format(parameters)


def test_heatmap_columns() -> None:
    """
    Test the heatmap columns plot method.
    """
    def _test_heatmap_columns_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.heatmap_columns(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [{'segment': [20, 40], 'show': False, 'close': True}, {
        'segment': [30, 50], 'mask_selfsubstitutions': True, 'show': False, 'close': True
    }]

    # Assert
    for parameters in list_params:
        assert _test_heatmap_columns_output(
            parameters
        ) is False, "plot_heatmap_columns failed with {} parameters".format(parameters)


def test_plot_miniheatmap() -> None:
    """
    Test the miniheatmap method.
    """

    # Define aux function
    def _test_plot_miniheatmap(parameters: dict) -> bool:
        error = False
        try:
            OBJ_TEST.miniheatmap(offset=0, **parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            error = True
        return error

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [
        {
            'show': False,
            'close': True,
        },
        {
            'mask_selfsubstitutions': True,
            'neworder_aminoacids': list('ACDEFGHIKLMNPQRSTVWY*'),
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
            'close': True,
        },
    ]

    # Assert
    for parameters in list_params:
        assert _test_plot_miniheatmap( # Assert that that set of parameters works on that object
            parameters,
        ) is False, "plot_miniheatmap failed with {} parameters".format(
            parameters,
        )
