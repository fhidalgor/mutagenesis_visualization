#!/usr/bin/env python
# coding: utf-8

# # Create objects

# This notebook will create objects with the datasets that are stored in the folder data. Having this objects easily accessible may pose useful for testing purposes.

# ## Import

# In[ ]:


import numpy as np
import pandas as pd
import os

# Local imports
try:
    from mutagenesis_visualization.main.scripts.code_class import Screen
    import mutagenesis_visualization.main.scripts.code_utils as code_utils
except ModuleNotFoundError:
    import import_notebook
    from code_class import Screen
    import code_utils


# ## HRas RBD

# In[ ]:


def hras_RBD():
    """
    Create object hras_RBD.
    """    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'HRas166_RBD.csv')
    except NameError:
        my_file = os.path.join('../../data', 'HRas166_RBD.csv')

    # Load enrichment scores
    hras_enrichment_RBD = np.genfromtxt(my_file, delimiter=',')

    # Define protein sequence
    hras_sequence = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

    # Define secondary structure
    secondary = [['L0'], ['β1'] * (9 - 1), ['L1'] * (15 - 9),
                 ['α1'] * (25 - 15), ['L2'] * (36 - 25), ['β2'] * (46 - 36),
                 ['L3'] * (48 - 46), ['β3'] * (58 - 48), ['L4'] * (64 - 58),
                 ['α2'] * (74 - 64), ['L5'] * (76 - 74), ['β4'] * (83 - 76),
                 ['L6'] * (86 - 83), ['α3'] * (103 - 86), ['L7'] * (110 - 103),
                 ['β5'] * (116 - 110), ['L8'] * (126 - 116),
                 ['α4'] * (137 - 126), ['L9'] * (140 - 137),
                 ['β6'] * (143 - 140), ['L10'] * (151 - 143),
                 ['α5'] * (172 - 151), ['L11'] * (190 - 172)]

    # Create object
    hras_RBD = Screen(
        dataset=hras_enrichment_RBD,
        sequence=hras_sequence,
        secondary=secondary
    )
    return hras_RBD


# ## Beta Lactamase

# In[ ]:


def bla_obj():
    """
    Create object for the beta lactamase dataset.
    """

    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_bla_raw.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_bla_raw.pkl')

    # Load data
    df_bla_raw = pd.read_pickle(my_file)

    # Parse
    df_bla, sequence_bla = code_utils.parse_pivot(
        df_bla_raw, col_data='DMS_amp_625_(b)'
    )

    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list(df_bla.index)
    neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW')

    # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
    start_position = df_bla.columns[0]

    # Define sequence. If you dont know the start of the sequence, just add X's
    sequence_bla_x = 'MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRP'        + 'EERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVREL'        + 'CSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTM'        + 'PAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGS'        + 'RGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW'

    # Define secondary structure
    secondary_bla = [['L0'] * 23, ['α1'] * (38 - 23), ['L1'] * 2,
                     ['β1'] * (48 - 40), ['L2'] * 5, ['β2'] * (57 - 53),
                     ['L3'] * (68 - 57), ['α2'] * (84 - 68), ['L4'] * (95 - 84),
                     ['α3'] * (100 - 95), ['L5'] * (103 - 100), ['α4'] *
                     (110 - 103), ['L6'] * (116 - 110), ['α5'] * (140 - 116),
                     ['L7'] * (1), ['α6'] * (153 - 141), ['L8'] * (164 - 153),
                     ['α7'] * (169 - 164), ['L9'] * (179 - 169),
                     ['α8'] * (194 - 179), ['L10'] * 3, ['α9'] * (210 - 197),
                     ['L11'] * (227 - 210), ['β3'] * (235 - 227),
                     ['L12'] * (240 - 235), ['β4'] * (249 - 240),
                     ['L13'] * (254 - 249), ['β5'] * (262 - 254),
                     ['L14'] * (266 - 262), ['α10'] * (286 - 266)]

    # Create objects
    bla_obj = Screen(
        df_bla, sequence_bla_x, aminoacids, start_position, 0, secondary_bla
    )
    return bla_obj


# ## Sumo1

# In[ ]:


def sumo_obj():
    """
    Create object for the sumo1 dataset.
    """
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_sumo1_raw.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_sumo1_raw.pkl')

    df_sumo1_raw = pd.read_pickle(my_file)

    # Parse
    df_sumo1, sequence_sumo1 = code_utils.parse_pivot(
        df_sumo1_raw, col_data='DMS'
    )

    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list(df_sumo1.index)
    neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW')

    # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
    start_position = df_sumo1.columns[0]

    # Full sequence
    sequence_sumo1 = 'MSDQEAKPSTEDLGDKKEGEYIKLKVIGQDSSEIHFKVKMTTHLKKLKESYCQRQGVPMN'        + 'SLRFLFEGQRIADNHTPKELGMEEEDVIEVYQEQTGGHSTV'
    # Define secondary structure
    secondary_sumo1 = [['L0'] * (20), ['β1'] * (28 - 20), ['L1'] * 3,
                       ['β2'] * (39 - 31), ['L2'] * 4, ['α1'] * (55 - 43),
                       ['L3'] * (6), ['β3'] * (65 - 61), ['L4'] * (75 - 65),
                       ['α2'] * (80 - 75), ['L5'] * (85 - 80),
                       ['β4'] * (92 - 85), ['L6'] * (101 - 92)]

    # Create objects
    sumo_obj = Screen(
        df_sumo1, sequence_sumo1, aminoacids, start_position, 1, secondary_sumo1
    )
    return sumo_obj


# ## MAPK1

# In[ ]:


def mapk1_obj():
    """
    Create object for the mapk1 dataset.
    """
    
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_mapk1_raw.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_mapk1_raw.pkl')

    # Read excel file
    df_mapk1_raw = pd.read_pickle(my_file)

    # Parse
    df_mapk1, sequence_mapk1 = code_utils.parse_pivot(
        df_mapk1_raw, col_data='DMS_DOX'
    )

    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list(df_mapk1.index)
    neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW')

    # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
    start_position = df_mapk1.columns[0]

    # Full sequence
    sequence_mapk1_x = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVRVAIK'                    +'KISPFEHQTYCQRTLREIKILLRFRHENIIGINDIIRAPTIEQMKDVYIVQDLMETDLYKLLKTQ'                    +'HLSNDHICYFLYQILRGLKYIHSANVLHRDLKPSNLLLNTTCDLKICDFGLARVADPDHDHTGFL'                    +'TEYVATRWYRAPEIMLNSKGYTKSIDIWSVGCILAEMLSNRPIFPGKHYLDQLNHILGILGSPSQ'                    +'EDLNCIINLKARNYLLSLPHKNKVPWNRLFPNADSKALDLLDKMLTFNPHKRIEVEQALAHPYLE'                    +'QYYDPSDEPIAEAPFKFDMELDDLPKEKLKELIFEETARFQPGYRS'

    # Create objects
    mapk1_obj = Screen(
        df_mapk1, sequence_mapk1_x, aminoacids, start_position, 0
    )
    return mapk1_obj


# ## UBE2I

# In[ ]:


def ube2i_obj():
    """
    Create object for the ube2i dataset.
    """
    
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_ube2i_raw.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_ube2i_raw.pkl')

    # Read excel file
    col_data = 'DMS'
    df_ube2i_raw = pd.read_pickle(my_file)

    # Parse
    df_ube2i, sequence_ube2i = code_utils.parse_pivot(
        df_ube2i_raw, col_data=col_data
    )

    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list(df_ube2i.index)
    neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW')

    # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
    start_position = df_ube2i.columns[0]  # Create object2i.columns[0]

    # Full sequence
    sequence_ube2i_x = 'MSGIALSRLAQERKAWRKDHPFGFVAVPTKNPDGTMNLMNWECAIPGKKGTP'                        +'WEGGLFKLRMLFKDDYPSSPPKCKFEPPLFHPNVYPSGTVCLSILEEDKDWRPAITIKQ'                        +'ILLGIQELLNEPNIQDPAQAEAYTIYCQNRVEYEKRVRAQAKKFAPS'

    # Define secondary structure
    secondary_ube2i = [
        ['α1'] * (20 - 1), ['L1'] * (24 - 20), ['β1'] * (30 - 24), ['L2'] * 5,
        ['β2'] * (46 - 35), ['L3'] * (56 - 46), ['β3'] * (63 - 56),
        ['L4'] * (73 - 63), ['β4'] * (77 - 73), ['L5'] * (93 - 77),
        ['α2'] * (98 - 93), ['L6'] * (107 - 98), ['α3'] * (122 - 107),
        ['L7'] * (129 - 122), ['α4'] * (155 - 129), ['L8'] * (160 - 155)
    ]

    # Create objects
    ube2i_obj = Screen(
        df_ube2i, sequence_ube2i_x, aminoacids, start_position, 1,
        secondary_ube2i
    )

    return ube2i_obj


# ## TAT

# In[ ]:


def tat_obj():
    """
    Create object for the tat dataset.
    """
    
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_tat.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_tat.pkl')

    # Read data
    df_tat = pd.read_pickle(my_file)

    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list(df_tat.index)
    neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW*')

    # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
    start_position = df_tat.columns[0]

    # Full sequence
    sequence_tat = 'MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRAHQ'                        +'NSQTHQASLSKQPTSQPRGDPTGPKE'

    # Define secondary structure
    secondary_tat = [['L1'] * (8), ['α1'] * (13 - 8), ['L2'] * (28 - 14),
                     ['α2'] * (41 - 28), ['L3'] * (90 - 41)]

    tat_obj = Screen(
        df_tat, sequence_tat, aminoacids, start_position, 0, secondary_tat
    )

    return tat_obj


# ## REV

# In[ ]:


def rev_obj():
    """
    Create object for the rev dataset.
    """
    
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_rev.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_rev.pkl')

    # Read dataset
    df_rev = pd.read_pickle(my_file)

    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list(df_rev.index)
    neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW*')

    # First residue of the hras_enrichment dataset. Because 1-Met was not mureved, the dataset starts at residue 2
    start_position = df_rev.columns[0]

    # Full sequence
    sequence_rev = 'MAGRSGDSDEDLLKAVRLIKFLYQSNPPPNPEGTRQARRNRRRRWRERQRQIHSISERIL'                    + 'STYLGRSAEPVPLQLPPLERLTLDCNEDCGTSGTQGVGSPQILVESPTILESGAKE'

    # Define secondary structure
    secondary_rev = [['L1'] * (8), ['α1'] * (25 - 8), ['L2'] * (33 - 25),
                     ['α2'] * (68 - 33), ['L3'] * (116 - 41)]

    rev_obj = Screen(
        df_rev, sequence_rev, aminoacids, start_position, 0, secondary_rev
    )
    return rev_obj


# ## α-synuclein

# In[ ]:


def asynuclein_obj():
    """
    Create object for the synuclein dataset.
    """
    
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_asynuclein.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_asynuclein.pkl')
    # Read dataset
    df_asynuclein = pd.read_pickle(my_file)

    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list(df_asynuclein.index)
    neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW')

    # First residue of the hras_enrichment dataset. Because 1-Met was not mureved, the dataset starts at residue 2
    start_position = df_asynuclein.columns[0]

    # Full sequence
    sequence_asynuclein = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK'                    + 'EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP'                    + 'DNEAYEMPSEEGYQDYEPEA'

    # Define secondary structure
    secondary_asynuclein = [['L1'] * (1), ['α1'] * (37 - 1), ['L2'] * (44 - 37),
                            ['α2'] * (92 - 44), ['L3'] * (140 - 92)]

    asynuclein_obj = Screen(
        df_asynuclein, sequence_asynuclein, aminoacids, start_position, 0,
        secondary_asynuclein
    )

    return asynuclein_obj


# ## APH(3) II

# In[ ]:


def aph_obj():
    """
    Create object for the aph dataset.
    """
    
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_aph.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_aph.pkl')
    # Read dataset
    df_aph = pd.read_pickle(my_file)

    aminoacids = list(df_aph.index)
    neworder_aminoacids = list('DEKHRGNQASTPCVYMILFW')

    # First residue of the hras_enrichment dataset. Because 1-Met was not mureved, the dataset starts at residue 2
    start_position = df_aph.columns[0]

    # Full sequence
    sequence_aph = 'MIEQDGLHAGSPAAWVERLFGYDWAQQTIGCSDAAVFRLSAQGRPVLFVKTDLSGALNELQ'                    + 'DEAARLSWLATTGVPCAAVLDVVTEAGRDWLLLGEVPGQDLLSSHLAPAEKVSIMADAMRR'                    + 'LHTLDPATCPFDHQAKHRIERARTRMEAGLVDQDDLDEEHQGLAPAELFARLKARMPDGED'                    + 'LVVTHGDACLPNIMVENGRFSGFIDCGRLGVADRYQDIALATRDIAEELGGEWADRFLVLY'                    + 'GIAAPDSQRIAFYRLLDEFF'

    # Define secondary structure
    secondary_aph = [['L1'] * (11), ['α1'] * (16 - 11), ['L2'] * (22 - 16),
                     ['β1'] * (26 - 22), ['L3'] * (34 - 26), ['β2'] * (40 - 34),
                     ['L4'] * (46 - 40), ['β3'] * (52 - 46), ['L5'] * (58 - 52),
                     ['α2'] * (72 - 58), ['L6'] * (79 - 72), ['β4'] * (85 - 79),
                     ['L7'] * (89 - 85), ['β5'] * (95 - 89), ['L8'] * (99 - 95),
                     ['β6'] * (101 - 99), ['L9'] * (107 - 101), ['α3'] *
                     (131 - 107), ['L10'] * (135 - 131), ['α4'] * (150 - 135),
                     ['L11'] * (158 - 150), ['α5'] * (163 - 158),
                     ['L12'] * (165 - 163), ['α6'] * (177 - 165),
                     ['L13'] * (183 - 177), ['β7'] * (187 - 183),
                     ['L14'] * (191 - 187), ['α7'] * (194 - 191), ['L15'] * (1),
                     ['β8'] * (199 - 195), ['L16'] * (201 - 199),
                     ['β9'] * (206 - 201), ['L17'] * (212 - 206),
                     ['β10'] * (216 - 212), ['α8'] * (245 - 216), ['L18'] * (4),
                     ['α9'] * (264 - 249)]

    aph_obj = Screen(
        np.log10(df_aph), sequence_aph, aminoacids, start_position, 0,
        secondary_aph
    )

    return aph_obj


# ## b11L5F

# In[ ]:


def b11L5F_obj():
    """
    Create object for the aph dataset.
    """
    
    # Use relative file import to access the data folder
    try:
        location = os.path.dirname(os.path.realpath(__file__))
        my_file = os.path.join(location, '../../data', 'df_b11L5F_raw.pkl')
    except NameError:
        my_file = os.path.join('../../data', 'df_b11L5F_raw.pkl')
    # Read dataset
    df_b11L5F_raw = pd.read_pickle(my_file)

    # Parse
    col_data = 'relative_tryp_stability_score'
    df_b11L5F, sequence_b11L5F = code_utils.parse_pivot(
        df_b11L5F_raw, col_data=col_data
    )

    # Order of amino acid substitutions in the hras_enrichment dataset
    aminoacids = list(df_b11L5F.index)
    neworder_aminoacids = list('DEKHRGNQASTPVYMILFW')

    # Sequence
    sequence_b11L5F = 'CRAASLLPGTWQVTMTNEDGQTSQGQMHFQPRSPYTLDVKAQGTISDGRPI'                        +'SGKGKVTCKTPDTMDVDITYPSLGNMKVQGQVTLDSPTQFKFDVTTSDGSKVTGTLQRQE'

    # First residue of the hras_enrichment dataset. Because 1-Met was not mureved, the dataset starts at residue 2
    start_position = df_b11L5F.columns[0]

    b11L5F_obj = Screen(
        df_b11L5F, sequence_b11L5F, aminoacids, start_position, 0
    )

    return b11L5F_obj

