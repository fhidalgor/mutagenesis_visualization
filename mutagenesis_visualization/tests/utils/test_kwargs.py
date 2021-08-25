"""
This module tests that kwargs work.
"""
from matplotlib.colors import LinearSegmentedColormap
from mutagenesis_visualization.main.utils.kwargs import generate_default_kwargs, _generate_colormap


#test that kwargs outputs a default dictionary
def test_kwargs() -> None:
    """
    Test the default kwargs generator.
    """
    assert isinstance(generate_default_kwargs(), dict), 'Error the kwarg type is not a dictionary.'


#test that _generatecolormap creates a matplotlib color map
def test_generatecolormap() -> None:
    """
    Test the aux function to generate a default color map used in heatmaps.
    """
    assert isinstance(
        _generate_colormap(), LinearSegmentedColormap
    ), 'Error the colormap type is not a matplotlib color map.'
