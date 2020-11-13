import code_utils

def code_utils._is_DNA():
    df = pd.DataFrame([[1,7,6,2],[0,0,0,0],[10,10,10,10],[1,1,1,1]])
    result = code_heatmaps._hierarchical_sort(df.T)
    assert (result == [2,0,1,3]).all()

