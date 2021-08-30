"""
This module contains the plotly figures tests.
"""

import traceback
import logging
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects
from mutagenesis_visualization.main.utils.data_paths import PDB_5P21

DEMO_OBJECTS: DemoObjects = DemoObjects()
OBJ_TEST = DEMO_OBJECTS.hras_rbd


def test_plotly_heatmap() -> None:
    """
    Test test_plotly_heatmap.
    """
    def _test_plotly_heatmap_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.plotly_heatmap(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except
            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: list = [
        {'show': False},
        {
            'show': False,
            'close': True,
            'title': 'Changing this parameter for fun',
        },
    ]

    # Assert
    for parameters in list_params:
        assert _test_plotly_heatmap_output(
            parameters
        ) is False, "plotly_heatmap failed with {} parameters".format(parameters)


def test_plotly_scatter() -> None:
    """
    Test plotly scatter plot.
    """
    def _test_plotly_scatter_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.plotly_scatter(
                OBJ_TEST, **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: list = [
        {'show': False},
        {
            'show': False,
            'close': True,
            'title': 'Changing this parameter for fun',
        },
        {'mode': 'mean', 'show': False, 'close': True, 'title': 'hello world'},
        {'mode': 'pointmutant', 'show': False, 'close': True, 'title': 'go bears'},
    ]

    # Assert
    for parameters in list_params:
        assert _test_plotly_scatter_output(
            parameters
        ) == False, "plotly_scatter failed with {} parameters".format(parameters)


def test_plotly_rank() -> None:
    """
    Test plotly rank plot.
    """
    def _test_plotly_rank_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.plotly_rank(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: list = [
        {'show': False},
        {
            'show': False,
            'close': True,
            'title': 'Changing this parameter for fun',
        },
        {'mode': 'mean', 'show': False, 'close': True, 'title': 'hello world'},
        {'mode': 'pointmutant', 'show': False, 'close': True, 'title': 'go bears'},
    ]

    # Assert
    for parameters in list_params:
        assert _test_plotly_rank_output(
            parameters
        ) == False, "plotly_rank failed with {} parameters".format(parameters)


def test_plotly_histogram() -> None:
    """
    Test plotly histogram plot.
    """
    def _test_plotly_histogram_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.plotly_histogram(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: list = [
        {'show': False},
        {
            'show': False,
            'close': True,
            'title': 'Changing this parameter for fun',
        },
        {'mode': 'mean', 'show': False, 'close': True, 'title': 'hello world'},
        {'mode': 'pointmutant', 'show': False, 'close': True, 'title': 'go bears'},
    ]

    # Assert
    for parameters in list_params:
        assert _test_plotly_histogram_output(
            parameters
        ) == False, "plotly_histogram failed with {} parameters".format(parameters)


def test_plotly_enrichment_bar() -> None:
    """
    Test plotly mean plot.
    """
    def _test_plotly_enrichment_bar_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.plotly_enrichment_bar(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: list = [
        {'show': False},
        {
            'show': False,
            'close': True,
            'title': 'Changing this parameter for fun',
        },
        {'mode': 'mean', 'show': False, 'close': True, 'title': 'hello world'},
        {'mode': 'A', 'show': False, 'close': True, 'title': 'go bears'},
    ]

    # Assert
    for parameters in list_params:
        assert _test_plotly_enrichment_bar_output(
            parameters
        ) == False, "plotly_enrichment_bar failed with {} parameters".format(parameters)


def test_plotly_scatter_3d() -> None:
    """
    Test plotly scater 3d.
    """
    def _test_plotly_scatter_3d_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.plotly_scatter_3d(**parameters)  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: list = [
        {
            'mode': 'A',
            'pdb_path': PDB_5P21,
            'title': 'Scatter 3d',
            'squared': True,
            'x_label': 'x',
            'y_label': 'y',
            'z_label': 'z',
            'show': False,
            'close': True,
        },
        {
            'mode': 'mean',
            'pdb_path': PDB_5P21,
            'title': 'Scatter 3d',
            'squared': False,
            'x_label': 'x',
            'y_label': 'y',
            'z_label': 'z',
            'show': False,
            'close': True,
        },
    ]

    # Assert
    for parameters in list_params:
        assert _test_plotly_scatter_3d_output(
            parameters
        ) == False, "plotly_scatter_3d failed with {} parameters".format(parameters)


def test_plotly_scatter_3d_pdbprop() -> None:
    """
    Test plotly scater 3d pdb properties.
    """
    def _test_plotly_scatter_3d_pdbprop_output(parameters: dict) -> bool:
        try:
            OBJ_TEST.plotly_scatter_3d_pdbprop(
                **parameters
            )  # pass dictionary as arguments of method
        except Exception as e:  # pylint: disable=broad-except

            logging.exception(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params: list = [
        {
            'mode': 'A',
            'plot': ['Distance', 'SASA', 'log B-factor'],
            'pdb_path': PDB_5P21,
            'show': False,
            'close': True,
        },
        {
            'plot': ['Distance', 'SASA', 'log B-factor'],
            'pdb_path': PDB_5P21,
            'show': False,
            'close': True,
        },
    ]

    # Assert
    for parameters in list_params:
        assert _test_plotly_scatter_3d_pdbprop_output(
            parameters
        ) == False, "plotly_scatter_3d_pdbprop failed with {} parameters".format(parameters)
