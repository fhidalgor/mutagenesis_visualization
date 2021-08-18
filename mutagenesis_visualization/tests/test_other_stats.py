"""
This module tests other stats folder.
"""
import traceback
import pandas as pd
import numpy as np
from pandas.core.frame import DataFrame

from mutagenesis_visualization.main.demo.demo_objects import DemoObjects

DEMO_OBJECTS:DemoObjects = DemoObjects()
OBJ_TEST_1 = DEMO_OBJECTS.hras_rbd
OBJ_TEST_2 = DEMO_OBJECTS.bla
DICT_OBJ: dict = {
    'obj_test_1': OBJ_TEST_1,
    'obj_test_2': OBJ_TEST_2
}

def test_rank():
    """
    Test the rank method.
    """

    # Define aux function
    def _test_rank(obj_test, parameters):
        try:
            obj_test.rank(
                **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {
            'show': False,
        },
        {
            'mode':'mean',
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in DICT_OBJ.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_rank( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "rank failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )

def test_roc():
    """
    Test the roc auc method.
    """
    # fake data for obj_test_1 (change for other objects)
    df_freq: DataFrame = pd.DataFrame()
    df_freq['Variant'] = OBJ_TEST_1.dataframe['Variant']
    df_freq['Class'] = np.random.randint(2, size=len(OBJ_TEST_1.dataframe))


    # Define aux function
    def _test_roc(obj_test, parameters):
        try:
            obj_test.roc(df_freq,
                **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {
            'show': False,
        },
        {
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in DICT_OBJ.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_roc( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "roc failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )


def test_cumulative():
    """
    Test the cumulative method.
    """

    # Define aux function
    def _test_cumulative(obj_test, parameters):
        try:
            obj_test.cumulative(
                **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:
            print(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {
            'show': False,
        },
        {
            'mode':'nonSNV',
            'show': False,
        },        {
            'mode':'SNV',
            'figsize': (3, 2.5),
            'y_label': r'$∆E^i_x$',
            'show': False,
        },
    ]

    # Assert
    for obj_label, obj_test in DICT_OBJ.items(): # Loop over the dictionary
        for parameters in list_params: # Loop over the parameters
            assert _test_cumulative( # Assert that that set of parameters works on that object
                obj_test,
                parameters,
            ) is False, "cumulative failed with {} object and {} parameters".format(
                obj_label,
                parameters,
            )
