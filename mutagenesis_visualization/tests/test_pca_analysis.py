"""
This module will test the pca analysis methods.
"""
import traceback
import logging
from typing import List
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects
from mutagenesis_visualization.main.classes.screen import Screen

DEMO_OBJECTS: DemoObjects = DemoObjects()
OBJ_TEST_1 = DEMO_OBJECTS.hras_rbd
OBJ_TEST_2 = DEMO_OBJECTS.sumo


def test_correlation() -> None:
    """
    Test the correlation method.
    """
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': OBJ_TEST_1,
        'obj_test_2': OBJ_TEST_2,
    }

    # Define aux function
    def _test_correlation(obj_test: Screen, parameters: dict) -> bool:
        try:
            obj_test.correlation(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except
            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [
        {
            'show': False,
            'close': True,
        },
        {
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
            'close': True,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items():  # Loop over the dictionary
        for parameters in list_params:  # Loop over the parameters
            assert _test_correlation( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "plot_correlation failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


def test_individual_correlation() -> None:
    """
    Test the individual correlation method.
    """
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': OBJ_TEST_1,
        'obj_test_2': OBJ_TEST_2,
    }

    # Define aux function
    def _test_individual_correlation(obj_test: Screen, parameters: dict) -> bool:
        try:
            obj_test.individual_correlation(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except
            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: List[dict] = [
        {
            'show': False,
            'close': True,
        },
        {
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
            'close': True,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items():  # Loop over the dictionary
        for parameters in list_params:  # Loop over the parameters
            assert _test_individual_correlation( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "plot_individual_correlation failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


def test_pca() -> None:
    """
    Test the pca method.
    """
    # Create dictionary with objects
    dict_obj = {
        'obj_test_1': OBJ_TEST_1,
        'obj_test_2': OBJ_TEST_2,
    }

    # Define aux function
    def _test_pca(obj_test: Screen, parameters: list) -> bool:
        try:
            obj_test.pca(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except
            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: list = [
        {
            'show': False,
            'close': True,
        },
        {
            'mode': 'individual',
            'show': False,
            'close': True,
        },
        {
            'mode': 'secondary',
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
            'close': True,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_obj.items():  # Loop over the dictionary
        for parameters in list_params:  # Loop over the parameters
            assert _test_pca( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "plot_pca failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )
