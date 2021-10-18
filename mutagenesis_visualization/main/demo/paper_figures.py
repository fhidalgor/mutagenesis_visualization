"""
In this module we will produce the figures that went into the paper.
"""
from typing import List, Tuple
import numpy as np
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects
from mutagenesis_visualization.main.classes.screen import Screen

HRAS_RBD = DemoObjects().hras_rbd
HRAS_COUNTS = DemoObjects().hras_counts
PATH: str = "docs/images/exported_images/paper/fig1{}.png"
ANC_RAS_ENRICHMENT = np.genfromtxt("figures_jk/enrichments/old/Ancestral167_RBD.csv", delimiter=',')
SHOW = False


def create_anc_ras_object() -> Tuple[Screen, List[Tuple[int, int]]]:
    """
    This dataset is not available since it is not published yet.
    """
    ancras: str = 'MTEYKLVVVGGGGVGKSALTIQFIQSHFVDEYDPTIEDSYRKQVVIDDEVAILDILDTAGQEEYSAMREQYMRNGEGFLLVYSITDRSSFDEISTYHEQILRVKDTDDVPMVLVGNKADLESRAVSMQEGQNLAKQLNVPFIETSAKQRMNVDEAFYTLVRVVRRH'
    hras: str = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQH'
    secondary: list = [['β1'] * (9 - 1), ['L1'] * (16 - 9), ['α1'] * (25 - 16), ['L2'] * (37 - 25),
                       ['β2'] * (46 - 37), ['L3'] * (49 - 46), ['β3'] * (58 - 49),
                       ['L4'] * (65 - 58), ['α2'] * (74 - 65), ['L5'] * (77 - 74),
                       ['β4'] * (83 - 77), ['L6'] * (87 - 83), ['α3'] * (103 - 87),
                       ['L7'] * (111 - 103), ['β5'] * (116 - 111), ['L8'] * (127 - 116),
                       ['α4'] * (137 - 127), ['L9'] * (141 - 137), ['β6'] * (143 - 141),
                       ['L10'] * (152 - 143), ['α5'] * (173 - 152), ['L11'] * (189 - 173)]
    kwargs: dict = {'secondary': secondary}
    aminoacids: List[str] = list('ACDEFGHIKLMNPQRSTVWY*')

    map_sequence_changes: List[Tuple[int, int]] = []
    for index, (i, j) in enumerate(zip(hras, ancras), start=1):
        if not i == j:
            map_sequence_changes.append((index, index))
    return Screen(ANC_RAS_ENRICHMENT, ancras, aminoacids, **kwargs), map_sequence_changes


ANC_RAS, MAP_SEQUENCE_CHANGES = create_anc_ras_object()


def generate_figures() -> None:
    """
    This figure will generate the figures that go into the paper.
    """
    HRAS_RBD.heatmap(title="", mask_selfsubstitutions=True,show_cartoon=True, output_file=PATH.format("a"), show=SHOW, dpi = 1200)
    HRAS_RBD.miniheatmap(mask_selfsubstitutions=True,title="Wt residue", output_file=PATH.format("b"), show=SHOW, dpi = 1200)
    HRAS_RBD.correlation(
        title="Amino acid correlation",
        output_file=PATH.format("c"),
        colorbar_scale=(0.5, 1),
        show=SHOW, dpi = 1200
    )
    HRAS_COUNTS.library_representation(title="", output_file=PATH.format("d"), show=SHOW, dpi = 1200)
    HRAS_RBD.enrichment_bar(
        title="", show_cartoon=True, output_file=PATH.format("e"), figsize=(6, 2.5), show=SHOW, dpi = 1200
    )
    HRAS_RBD.pca(
        adjust_labels=True, title="", mode="secondary", output_file=PATH.format("f"), show=SHOW, dpi = 1200
    )
    HRAS_RBD.sequence_differences(
        ANC_RAS,
        MAP_SEQUENCE_CHANGES,
        legend_labels=("H-Ras to Anc-Ras", "Anc-Ras to H-Ras"),
        bins=10,
        title="",
        output_file=PATH.format("h"),
        show=SHOW, dpi = 1200
    )

    # Create noise for class
    mu, sigma = 0.2, 0.5
    df_class = HRAS_RBD.dataframes.df_notstopcodons[-1].copy()
    np.random.seed(1)
    noise = np.random.normal(mu, sigma, len(df_class["Score"]))
    df_class["Score+Noise"] = df_class["Score"] + noise
    df_class["Class"] = df_class["Score+Noise"].astype(int)
    df_class.loc[df_class["Class"] < 0] = 0
    df_class.loc[df_class["Class"] > 1] = 1
    HRAS_RBD.roc(df_class, mode="mean", title="", output_file=PATH.format("g"), show=SHOW, dpi = 1200)


if __name__ == "__main__":
    generate_figures()
