from typing import List
from mutagenesis_visualization.main.classes.counts import Counts
from mutagenesis_visualization.main.classes.screen import Screen
from mutagenesis_visualization.main.demo.demo_data import load_demo_datasets, return_hras_counts
from mutagenesis_visualization.main.demo.demo_figures import run_demo
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects


def test_run_demo():
    """
    This function will test that demo is capable of generating the
    types of figures ('heatmap', 'miniheatmap', 'mean', 'kernel', 'pca',
    'position', 'secondary_mean', 'correlation', 'individual_correlation')
    that demo() is supposed to. Will raise an error if at least one of
    the plots does not work.
    """
    def _test_output_run_demo(argument: str) -> None:
        """
        Aux function for test_demo.
        Will try to run a demo function, will return True if there is an error.
        """
        try:
            run_demo(argument, show=False)
            return False
        except:
            return True

    arguments: List[str] = [
        'heatmap', 'miniheatmap', 'mean', 'kernel', 'pca', 'position',
        'secondary_mean', 'correlation', 'individual_correlation'
    ]
    solutions = [_test_output_run_demo(argument) for argument in arguments]

    failed = [arg for arg, sol in zip(arguments, solutions) if sol is True]

    assert any(solutions) == False, 'the following failed: {}'.format(
        str(failed)
    )

def test_load_demo_datasets() -> None:
    """
    Test that function returns dictionary.
    """
    assert (
        str(type(load_demo_datasets()))
    ) == "<class 'dict'>", "function demo_datasets failed"

def test_return_hras_counts() -> None:
    """
    Test return_hras_counts.
    """
    demo_object: Counts = return_hras_counts()
    assert isinstance(demo_object, Counts)

def test_demo_objects() -> None:
    """
    Test class DemoObjects.
    """
    demo_object: DemoObjects = DemoObjects()
    assert isinstance(demo_object.hras_RBD, Screen)
    assert isinstance(demo_object.bla, Screen)
    assert isinstance(demo_object.sumo, Screen)
    assert isinstance(demo_object.mapk1, Screen)
    assert isinstance(demo_object.ube2i, Screen)
    assert isinstance(demo_object.tat, Screen)
    assert isinstance(demo_object.rev, Screen)
    assert isinstance(demo_object.asynuclein, Screen)
    assert isinstance(demo_object.aph, Screen)
    assert isinstance(demo_object.b11L5F, Screen)
