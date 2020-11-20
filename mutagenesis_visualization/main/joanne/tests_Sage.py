Instructions for test_plotly

functions to test
	rank_plotly
	rank_plotly(mode='mean')
	scatter_plotly
	scatter_plotly(mode='mean')
	
	you can use this function copied from code_demo to generate an example (see below). Change the name to something like _test_plotly_generate_example and no parameters. Docstring would be different too. THen use the object to test if it produces the plots. Because it is difficult to test if a 
    plot was made, the only thing we should test is whether the function is able to work properly. See notebook test_demo.ipynb for an example. I copied the function that is used to test if a function runs properly (_test_output_demo). You will also have to do some changes.
	
def demo(figure='heatmap', show=True):
    """
    Performs a demonstration of the mutagenesis_visualization software.

    Parameters
    -----------
    figure : str, default 'heatmap'
        There are a few example plots that can be displayed to test the package is working on your station.
        The options are 'heatmap', 'miniheatmap', 'mean', 'kernel', 'pca'
        'position', 'secondary_mean', 'correlation', 'individual_correlation'. Check the documentation for more information.
    
    show : boolean, default True
        If True, will do plt.show() for each figure
        
    Returns
    -------
    None.
    """
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../data', 'HRas166_RBD.csv')
    except NameError:
        my_file = os.path.join('../data', 'HRas166_RBD.csv')

    # Load enrichment scores
    hras_enrichment_RBD = np.genfromtxt(my_file, delimiter=',')

    # Define protein sequence
    hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

    # Define secondary structure
    secondary = [['L0'], ['β1']*(9-1), ['L1']*(15-9), ['α1']*(25-15), ['L2']*(36-25), ['β2']*(46-36), ['L3']*(48-46),
                 ['β3']*(58-48), ['L4'] * (64-58), ['α2'] *
                 (74-64), ['L5']*(76-74), ['β4']*(83-76),
                 ['L6']*(86-83), ['α3']*(103-86), ['L7']*(110-103), ['β5'] *
                 (116-110), ['L8']*(126-116), ['α4']*(137-126),
                 ['L9']*(140-137), ['β6']*(143-140), ['L10']*(151-143), ['α5']*(172-151), ['L11']*(190-172)]

    # Create object
    hras_RBD = Screen(dataset=hras_enrichment_RBD,
                      sequence=hras_sequence, secondary=secondary)
    return hras_RBD

def _test_output_demo(argument):
    '''
    Aux function for test_demo.
    Will try to run a demo function, will return True if there is an error.
    '''
    error = False
    try:
        demo(argument, show=False)
    except:
        error = True
    return error
    
def test_demo():
    '''
    This function will test that demo is capable of generating the 
    types of figures ('heatmap', 'miniheatmap', 'mean', 'kernel', 'pca',
    'position', 'secondary_mean', 'correlation', 'individual_correlation') that demo()
    is supposed to. Will raise an error if at least one of the plots does not work.
    '''
    
    arguments = ['heatmap', 'miniheatmap', 'mean', 'kernel', 'pca',
                'position', 'secondary_mean', 'correlation', 'individual_correlation']
    solutions = [_test_output_demo(argument) for argument in arguments]
    assert any(solutions)==False, 'error when running the demo figures'