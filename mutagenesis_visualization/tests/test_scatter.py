"""
This module contains the scatter tests.
"""

import traceback
import logging
from typing import Dict
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects
from mutagenesis_visualization.main.classes.screen import Screen


def test_plot_scatter() -> None:
    """
    Test plot scatter function.
    """
    demo_objects: DemoObjects = DemoObjects()
    dict_objects: Dict[str, Screen] = {"hras": demo_objects.hras_rbd, "bla": demo_objects.bla}

    # Define aux function
    def _test_plot_scatter(obj_test: Screen, parameters: list) -> bool:
        try:
            obj_test.scatter(obj_test, **parameters)  # pass dictionary as arguments of method
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
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
            'close': True,
        },
    ]

    # Assert
    for obj_label, obj_test in dict_objects.items():  # Loop over the dictionary
        for parameters in list_params:  # Loop over the parameters
            assert _test_plot_scatter( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "plot_scatter failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )
