"""
This module tests the bar plots.
"""
import traceback
import logging
from typing import Dict, List
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects
from mutagenesis_visualization.main.classes.screen import Screen

DEMO_OBJECTS: DemoObjects = DemoObjects()
OBJ_TEST_1 = DEMO_OBJECTS.hras_rbd
OBJ_TEST_2 = DEMO_OBJECTS.bla
DICT_OBJ: Dict[str, Screen] = {'obj_test_1': OBJ_TEST_1, 'obj_test_2': OBJ_TEST_2}
OBJ_TEST_COUNTS = DEMO_OBJECTS.hras_counts


def test_enrichment_bar() -> None:
    """
    Test mean bar plot.
    """

    # Define aux function
    def _test_enrichment_bar_output(obj_test: Screen, parameters: dict) -> bool:

        try:
            obj_test.enrichment_bar(**parameters)  # pass dictionary as arguments of method
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
            'mode': 'A',
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
            'close': True,
            'return_object': True,
        },
    ]

    # Assert
    for obj_label, obj_test in DICT_OBJ.items():  # Loop over the dictionary
        for parameters in list_params:  # Loop over the parameters
            assert _test_enrichment_bar_output( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "enrichment_bar failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


def test_differential() -> None:
    """
    Test the differential bar plot.
    """
    def _test_differential_output(obj_test: Screen, parameters: dict) -> bool:
        try:
            obj_test.differential(obj_test, **parameters)  # pass dictionary as arguments of method
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
            'show_cartoon': True,
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
            'close': True,
        },
    ]

    # Assert
    for obj_label, obj_test in DICT_OBJ.items():  # Loop over the dictionary
        for parameters in list_params:  # Loop over the parameters
            assert _test_differential_output( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "differential failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


#does not work when the position is 1 --should make a note of this somewhere
def test_position_bar() -> None:
    """
    Test the position bar plot.
    """
    def _test_position_bar_output(obj_test: Screen, parameters: dict) -> bool:
        try:
            obj_test.position_bar(
                position=115, **parameters
            )  # pass dictionary as arguments of method
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
    for obj_label, obj_test in DICT_OBJ.items():  # Loop over the dictionary
        for parameters in list_params:  # Loop over the parameters
            assert _test_position_bar_output( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "position_bar failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


def test_secondary() -> None:
    """
    Test the secondary bar plot.
    """

    # Define aux function
    def _test_secondary(obj_test: Screen, parameters: dict) -> bool:
        try:
            obj_test.secondary_mean(**parameters)  # pass dictionary as arguments of method
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
    for obj_label, obj_test in DICT_OBJ.items():  # Loop over the dictionary
        for parameters in list_params:  # Loop over the parameters
            assert _test_secondary( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "secondary failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


def test_library_representation() -> None:
    """
    Test the library representation method.
    """
    def _test_library_representation(parameters: dict) -> bool:
        try:
            OBJ_TEST_COUNTS.library_representation(**parameters)
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
            'figsize': (
                3, 2.5
            ),  #legend does not change to scale --line 70-80 of code_bar of plot_library_representation
            'y_label':
            r'$∆E^i_x$',  #y_label does not change --line 47 of code_bar of plot_library_representation
            'show': False,
            'close': True,
            'title': 'go bears'
        },
    ]

    # Assert
    for parameters in list_params:  # Loop over the parameters
        assert _test_library_representation( # Assert that that set of parameters works on that object
            parameters,
        ) is False, "plot_library failed with {} parameters".format(
            parameters,
        )


def test_mean_counts() -> None:
    """
    Test the mean counts bar plot.
    """
    def _test_mean_counts_output(parameters: dict) -> bool:
        try:
            OBJ_TEST_COUNTS.mean_counts(**parameters)  # pass dictionary as arguments of method
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
            'figsize': (3, 2.5), 'y_label': r'$∆E^i_x$', 'show': False, 'close': True,
            'text_labels': [[1, 1, 'label']]
        },
    ]

    # Assert
    for parameters in list_params:  # Loop over the parameters
        assert _test_mean_counts_output( # Assert that that set of parameters works on that object
            parameters,
        ) is False, "plot_meancounts failed with {} parameters".format(
            parameters,
        )
