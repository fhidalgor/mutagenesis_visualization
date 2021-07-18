"""
This module tests the kernel and histogram methods.
"""
import traceback
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects

DEMO_OBJECTS:DemoObjects = DemoObjects()
OBJ_TEST = DEMO_OBJECTS.hras_rbd

def test_kernel():
    def _test_kernel_output(parameters):
        try:
            OBJ_TEST.kernel(
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
        {'show': False},
        {'cumulative': True, 'show': False},
        {'y_label': 'testing', 'title':'this is a test', 'x_label':'testing', 'xscale':2, 'show': False} #y_label does not change
    ]

    # Assert
    for parameters in list_params:
        assert _test_kernel_output(
            parameters
        ) == False, "plot_kernel failed with {} parameters".format(parameters)


def test_histogram():
    def _test_histogram_output(parameters):
        try:
            OBJ_TEST.histogram(
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
        {'show': False},
        {'population': 'SNV', 'show': False},
        {'population': 'nonSNV', 'show': False},
    ]

    # Assert
    for parameters in list_params:
        assert _test_histogram_output(
            parameters
        ) == False, "plot_histogram failed with {} parameters".format(
            parameters
        )

def test_multiple_kernel():
    """
    Test of the multiple kernels method.
    """
    def _test_multiple_kernel(parameters):
        try:
            OBJ_TEST.multiple_kernel(DEMO_OBJECTS.bla, **parameters)
        except Exception as e:

            print(e)
            print(traceback.format_exc())
            return True
        return False

    # Define dictionary of parameters
    # Each dict in the list will be a different set of parameters
    list_params = [
        {   'label_kernels': ["H-Ras", "BLA"],
            'show': False,
        },
        {   'label_kernels': ["H-Ras", "BLA"],
            'figsize': (3, 2.5), 'y_label': r'$∆E^i_x$', 'show': False, 'title':
            'go bears'
        },
    ]

    # Assert
    for parameters in list_params:  # Loop over the parameters
        assert _test_multiple_kernel( # Assert that that set of parameters works on that object
            parameters,
        ) is False, "plot_multiplekernel failed with {} object and {} parameters".format(
                parameters
            )
