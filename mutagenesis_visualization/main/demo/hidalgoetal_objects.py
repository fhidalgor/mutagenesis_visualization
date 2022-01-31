"""
In this module I will generate the objects from the Hidalgo et al.
2022 eLife paper.
"""
from typing import List, Dict
import numpy as np
from mutagenesis_visualization.main.utils.data_paths import DATA_FOLDER
from mutagenesis_visualization.main.classes.screen import Screen

# Load data
HRas166_GAP_r0 = np.genfromtxt(DATA_FOLDER + 'HRas166_GAP_r0.csv', delimiter=',')
HRas166_GAP_r1 = np.genfromtxt(DATA_FOLDER + 'HRas166_GAP_r1.csv', delimiter=',')

HRas166_GAPGEF_r0 = np.genfromtxt(DATA_FOLDER + 'HRas166_GAPGEF_r0.csv', delimiter=',')
HRas166_GAPGEF_r1 = np.genfromtxt(DATA_FOLDER + 'HRas166_GAPGEF_r1.csv', delimiter=',')

HRas166_RBD_r0 = np.genfromtxt(DATA_FOLDER + 'HRas166_RBD_r0.csv', delimiter=',')
HRas166_RBD_r1 = np.genfromtxt(DATA_FOLDER + 'HRas166_RBD_r1.csv', delimiter=',')
HRas166_RBD_r2 = np.genfromtxt(DATA_FOLDER + 'HRas166_RBD_r2.csv', delimiter=',')

HRas188_BaF3_r0 = np.genfromtxt(DATA_FOLDER + 'HRas188_BaF3_r0.csv', delimiter=',')
HRas188_BaF3_r1 = np.genfromtxt(DATA_FOLDER + 'HRas188_BaF3_r1.csv', delimiter=',')

HRas180_GAP_r0 = np.genfromtxt(DATA_FOLDER + 'HRas180_GAP_r0.csv', delimiter=',')
HRas180_GAP_r1 = np.genfromtxt(DATA_FOLDER + 'HRas180_GAP_r1.csv', delimiter=',')
HRas180_GAP_r2 = np.genfromtxt(DATA_FOLDER + 'HRas180_GAP_r2.csv', delimiter=',')
HRas180_GAP_r3 = np.genfromtxt(DATA_FOLDER + 'HRas180_GAP_r3.csv', delimiter=',')

HRas180_RBD_r0 = np.genfromtxt(DATA_FOLDER + 'HRas180_RBD_r0.csv', delimiter=',')
HRas180_RBD_r1 = np.genfromtxt(DATA_FOLDER + 'HRas180_RBD_r1.csv', delimiter=',')
HRas180_RBD_r2 = np.genfromtxt(DATA_FOLDER + 'HRas180_RBD_r2.csv', delimiter=',')
HRas180_RBD_r3 = np.genfromtxt(DATA_FOLDER + 'HRas180_RBD_r3.csv', delimiter=',')

KRas173_RBD_r0 = np.genfromtxt(DATA_FOLDER + 'KRas173_RBD_r0.csv', delimiter=',')
KRas173_RBD_r1 = np.genfromtxt(DATA_FOLDER + 'KRas173_RBD_r1.csv', delimiter=',')
KRas173_RBD_r2 = np.genfromtxt(DATA_FOLDER + 'KRas173_RBD_r2.csv', delimiter=',')
KRas173_RBD_r3 = np.genfromtxt(DATA_FOLDER + 'KRas173_RBD_r3.csv', delimiter=',')

KRas165_RBD_r0 = np.genfromtxt(DATA_FOLDER + 'KRas165_RBD_r0.csv', delimiter=',')
KRas165_RBD_r1 = np.genfromtxt(DATA_FOLDER + 'KRas165_RBD_r1.csv', delimiter=',')

KRas165_GAP_r0 = np.genfromtxt(DATA_FOLDER + 'KRas165_GAP_r0.csv', delimiter=',')
KRas165_GAP_r1 = np.genfromtxt(DATA_FOLDER + 'KRas165_GAP_r1.csv', delimiter=',')
KRas165_GAP_r2 = np.genfromtxt(DATA_FOLDER + 'KRas165_GAP_r2.csv', delimiter=',')

KRas173_GAP_r0 = np.genfromtxt(DATA_FOLDER + 'KRas173_GAP_r0.csv', delimiter=',')
KRas173_GAP_r1 = np.genfromtxt(DATA_FOLDER + 'KRas173_GAP_r1.csv', delimiter=',')
KRas173_GAP_r2 = np.genfromtxt(DATA_FOLDER + 'KRas173_GAP_r2.csv', delimiter=',')
KRas173_GAP_r3 = np.genfromtxt(DATA_FOLDER + 'KRas173_GAP_r3.csv', delimiter=',')

KRas165_GAPGEF_r0 = np.genfromtxt(DATA_FOLDER + 'KRas165_GAPGEF_r0.csv', delimiter=',')
KRas165_GAPGEF_r1 = np.genfromtxt(DATA_FOLDER + 'KRas165_GAPGEF_r1.csv', delimiter=',')

KRas173_GAPGEF_r0 = np.genfromtxt(DATA_FOLDER + 'KRas173_GAPGEF_r0.csv', delimiter=',')
KRas173_GAPGEF_r1 = np.genfromtxt(DATA_FOLDER + 'KRas173_GAPGEF_r1.csv', delimiter=',')

KRas165_GEF_r0 = np.genfromtxt(DATA_FOLDER + 'KRas165_GEF_r0.csv', delimiter=',')
KRas165_GEF_r1 = np.genfromtxt(DATA_FOLDER + 'KRas165_GEF_r1.csv', delimiter=',')

KRas173_GEF_r0 = np.genfromtxt(DATA_FOLDER + 'KRas173_GEF_r0.csv', delimiter=',')
KRas173_GEF_r1 = np.genfromtxt(DATA_FOLDER + 'KRas173_GEF_r1.csv', delimiter=',')

KRas173Q61L_RBD_r0 = np.genfromtxt(DATA_FOLDER + 'KRas173Q61L_RBD_r0.csv', delimiter=',')
KRas173Q61L_RBD_r1 = np.genfromtxt(DATA_FOLDER + 'KRas173Q61L_RBD_r1.csv', delimiter=',')
KRas173Q61L_RBD_r2 = np.genfromtxt(DATA_FOLDER + 'KRas173Q61L_RBD_r2.csv', delimiter=',')
KRas173Q61L_RBD_r3 = np.genfromtxt(DATA_FOLDER + 'KRas173Q61L_RBD_r3.csv', delimiter=',')

KRas173Q61L_GAP_r0 = np.genfromtxt(DATA_FOLDER + 'KRas173Q61L_GAP_r0.csv', delimiter=',')
KRas173Q61L_GAP_r1 = np.genfromtxt(DATA_FOLDER + 'KRas173Q61L_GAP_r1.csv', delimiter=',')
#KRas173Q61L_GAP_r2 = np.genfromtxt(DATA_FOLDER + 'KRas173Q61L_GAP_r2.csv', delimiter=',')

# Sequences
hras: str = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPGCMSCKCVLS'
kras: str = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKD'
ancras: str = 'MTEYKLVVVGGGGVGKSALTIQFIQSHFVDEYDPTIEDSYRKQVVIDDEVAILDILDTAGQEEYSAMREQYMRNGEGFLLVYSITDRSSFDEISTYHEQILRVKDTDDVPMVLVGNKADLESRAVSMQEGQNLAKQLNVPFIETSAKQRMNVDEAFYTLVRVVRRH'
srosettaras: str = 'MTEYRLVVVGTGGVGKSALTIQLIQQHFVTEYDPTIEDSYRKHVSIDDEACLLDILDTAGQEDYSAMRDQYMRTGEGFLCVYSIDSQQSLDEIHSFREQILRVKDQDEVPMILVGNKCDLEHREVSTEAGQAVAKSYSIPFMETSAKKRINVEEAFYQLVREIRKY'
krasQ61L: str = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGLEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHHYREQIKRVKDSEDVPMVLVGNKCDLPSRTVDTKQAQDLARSYGIPFIETSAKTRQGVDDAFYTLVREIRKHKEKMSKD'

secondary: list = [['β1'] * (9 - 1), ['L1'] * (16 - 9), ['α1'] * (25 - 16), ['L2'] * (37 - 25),
             ['β2'] * (46 - 37), ['L3'] * (49 - 46), ['β3'] * (58 - 49), ['L4'] * (65 - 58),
             ['α2'] * (74 - 65), ['L5'] * (77 - 74), ['β4'] * (83 - 77),
             ['L6'] * (87 - 83), ['α3'] * (103 - 87), ['L7'] * (111 - 103), ['β5'] * (116 - 111),
             ['L8'] * (127 - 116), ['α4'] * (137 - 127), ['L9'] * (141 - 137), ['β6'] * (143 - 141),
             ['L10'] * (152 - 143), ['α5'] * (173 - 152), ['L11']*(189-173) ]

kwargs: dict = {'secondary': secondary}
aminoacids: List[str] = list('ACDEFGHIKLMNPQRSTVWY*')

# Create objects
HRAS_166_GAP: Screen = Screen([HRas166_GAP_r0, HRas166_GAP_r1], hras, aminoacids, **kwargs)
HRAS_166_RBD: Screen = Screen([HRas166_RBD_r0, HRas166_RBD_r1, HRas166_RBD_r2], hras, aminoacids, **kwargs)
HRAS_166_GAPGEF: Screen = Screen([HRas166_GAPGEF_r0, HRas166_GAPGEF_r1], hras, aminoacids,  **kwargs)
HRAS_188_BAF3: Screen = Screen([HRas188_BaF3_r0, HRas188_BaF3_r1], hras, aminoacids,  **kwargs)
HRAS_180_GAP: Screen = Screen([HRas180_GAP_r0, HRas180_GAP_r1, HRas180_GAP_r3], hras, aminoacids,  **kwargs)
HRAS_180_RBD: Screen = Screen([HRas180_RBD_r0, HRas180_RBD_r1, HRas180_RBD_r2, HRas180_RBD_r3], hras, aminoacids,  **kwargs)
KRAS_165_GAP: Screen = Screen([KRas165_GAP_r0, KRas165_GAP_r1, KRas165_GAP_r2], kras, aminoacids,  **kwargs)
KRAS_165_RBD: Screen = Screen([KRas165_RBD_r0, KRas165_RBD_r1], kras, aminoacids,  **kwargs)
KRAS_173_GAP: Screen = Screen([KRas173_GAP_r0, KRas173_GAP_r1, KRas173_GAP_r2, KRas173_GAP_r3], kras, aminoacids,  **kwargs)
KRAS_173_RBD: Screen = Screen([KRas173_RBD_r0, KRas173_RBD_r1, KRas173_RBD_r2, KRas173_RBD_r3], kras, aminoacids,  **kwargs)
KRAS_165_GAPGEF: Screen = Screen([KRas165_GAPGEF_r0, KRas165_GAPGEF_r1], kras, aminoacids,  **kwargs)
KRAS_173_GAPGEF: Screen = Screen([KRas173_GAPGEF_r0, KRas173_GAPGEF_r1], kras, aminoacids,  **kwargs)
KRAS_165_GEF: Screen = Screen([KRas165_GEF_r0, KRas165_GEF_r1], kras, aminoacids,  **kwargs)
KRAS_173_GEF: Screen = Screen([KRas173_GEF_r0, KRas173_GEF_r1], kras, aminoacids,  **kwargs)
KRASQ61L_173_RBD: Screen = Screen([KRas173Q61L_RBD_r0, KRas173Q61L_RBD_r1, KRas173Q61L_RBD_r2, KRas173Q61L_RBD_r3], krasQ61L, aminoacids,  **kwargs)
KRASQ61L_173_GAP: Screen = Screen([KRas173Q61L_GAP_r0, KRas173Q61L_GAP_r1], krasQ61L, aminoacids,  **kwargs)

OBJECTS: Dict[str, Screen] = {
"H-Ras_2-166_GAP": HRAS_166_GAP,
"H-Ras_2-166_RBD": HRAS_166_RBD,
"H-Ras_2-188_BaF3": HRAS_188_BAF3,
"H-Ras_2-180_GAP": HRAS_180_GAP,
"H-Ras_2-180_RBD": HRAS_180_RBD,
"K-Ras_2-165_GAP": KRAS_165_GAP,
"K-Ras_2-165_GAPGEF": KRAS_165_GAPGEF,
"K-Ras_2-173_GAPGEF": KRAS_173_GAPGEF,
"K-Ras_2-165_GEF": KRAS_165_GEF,
"K-Ras_2-173_GEF": KRAS_173_GEF,
"K-Ras_2-165_RBD": KRAS_165_RBD,
"K-Ras_2-173_GAP": KRAS_173_GAP,
"K-Ras_2-173_RBD": KRAS_173_RBD,
"H-Ras_2-166_GAPGEF": HRAS_166_GAPGEF,
"K-Ras_2-173_Q61L_GAP":KRASQ61L_173_GAP,
"K-Ras_2-173_Q61L_RBD":KRASQ61L_173_RBD,
}

def prettify_title(old_title: str) -> str:
    """
    Will make the title beautiful
    """
    new_title = old_title.replace("_", " ")
    replacements = {" 2-165":"$^{2-165}$"," 2-166":"$^{2-166}$", " 2-173":"$^{2-173}$", " 2-180":"$^{2-180}$", " 2-188":"$^{1-188}$",
                    " GAPGEF":"+GAP+GEF", " GAP":"+GAP", " GEF":"+GEF", "BaF3": "in Ba/F3", " RBD":"", "S.Rosetta":"S. Rosetta"}
    for old, new in replacements.items():
        new_title = new_title.replace(old, new)
    return new_title