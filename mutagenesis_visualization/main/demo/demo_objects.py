"""
This module will host the DemoObjects class.
"""
from typing import Dict, Union, List
from pandas.core.frame import DataFrame
import numpy as np

from mutagenesis_visualization.main.classes.screen import Screen
from mutagenesis_visualization.main.demo.demo_data import load_demo_datasets

DEMO_DATASETS: Dict[str, Union[np.array, DataFrame]] = load_demo_datasets()


class DemoObjects:
    """
    This class will load demo objects.
    """
    def __init__(self) -> None:
        self.aph: Screen = self._generate_aph_obj()
        self.asynuclein: Screen = self._generate_asynuclein_obj()
        self.b11L5F: Screen = self._generate_b11L5F_obj()
        self.bla: Screen = self._generate_bla_obj()
        self.hras_rbd: Screen = self._generate_hras_rbd_object()
        self.mapk1: Screen = self._generate_mapk1_obj()
        self.rev: Screen = self._generate_rev_obj()
        self.sumo: Screen = self._generate_sumo_obj()
        self.tat: Screen = self._generate_tat_obj()
        self.ube2i: Screen = self._generate_ube2i_obj()

    def _generate_hras_rbd_object(self) -> Screen:
        """
        Create object hras_RBD.
        """
        # Define protein sequence
        hras_sequence: str = 'MTEYKLVVVGAGGVGKSALTIQLIQNHFVDEYDPTIEDSYRKQVVIDGETCLLDILDTAGQEEYSAMRDQYMRTGEGFLCVFAINNTKSFEDIHQYREQIKRVKDSDDVPMVLVGNKCDLAARTVESRQAQDLARSYGIPYIETSAKTRQGVEDAFYTLVREIRQHKLRKLNPPDESGPG'

        # Define secondary structure
        secondary: list = [
            ['L0'], ['β1'] * (9 - 1), ['L1'] * (15 - 9), ['α1'] * (25 - 15), ['L2'] * (36 - 25),
            ['β2'] * (46 - 36), ['L3'] * (48 - 46), ['β3'] * (58 - 48), ['L4'] * (64 - 58),
            ['α2'] * (74 - 64), ['L5'] * (76 - 74), ['β4'] * (83 - 76), ['L6'] * (86 - 83),
            ['α3'] * (103 - 86), ['L7'] * (110 - 103), ['β5'] * (116 - 110), ['L8'] * (126 - 116),
            ['α4'] * (137 - 126), ['L9'] * (140 - 137), ['β6'] * (143 - 140), ['L10'] * (151 - 143),
            ['α5'] * (172 - 151), ['L11'] * (190 - 172)
        ]

        return Screen(
            dataset=DEMO_DATASETS['array_hras_RBD'], sequence=hras_sequence, aminoacids = list('ACDEFGHIKLMNPQRSTVWY*'), secondary=secondary
        )

    def _generate_bla_obj(self) -> Screen:
        """
        Create object for the beta lactamase dataset.
        """
        # Order of amino acid substitutions in the hras_enrichment dataset
        aminoacids: List[str] = list(DEMO_DATASETS['df_bla'].index)

        # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_bla'].columns[0]

        # Define sequence. If you dont know the start of the sequence, just add X's
        sequence_bla_x = 'MSIQHFRVALIPFFAAFCLPVFAHPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRP' + 'EERFPMMSTFKVLLCGAVLSRVDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVREL' + 'CSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTM' + 'PAAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGS' + 'RGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW'

        # Define secondary structure
        secondary_bla = [['L0'] * 23, ['α1'] * (38 - 23), ['L1'] * 2, ['β1'] * (48 - 40),
                         ['L2'] * 5, ['β2'] * (57 - 53), ['L3'] * (68 - 57), ['α2'] * (84 - 68),
                         ['L4'] * (95 - 84), ['α3'] * (100 - 95), ['L5'] * (103 - 100),
                         ['α4'] * (110 - 103), ['L6'] * (116 - 110), ['α5'] * (140 - 116),
                         ['L7'] * (1), ['α6'] * (153 - 141), ['L8'] * (164 - 153),
                         ['α7'] * (169 - 164), ['L9'] * (179 - 169), ['α8'] * (194 - 179), ['L10'] *
                         3, ['α9'] * (210 - 197), ['L11'] * (227 - 210), ['β3'] * (235 - 227),
                         ['L12'] * (240 - 235), ['β4'] * (249 - 240), ['L13'] * (254 - 249),
                         ['β5'] * (262 - 254), ['L14'] * (266 - 262), ['α10'] * (286 - 266)]

        return Screen(
            DEMO_DATASETS['df_bla'], sequence_bla_x, aminoacids, start_position, 0, secondary_bla
        )

    def _generate_sumo_obj(self) -> Screen:
        """
        Create object for the sumo1 dataset.
        """
        # Order of amino acid substitutions in the hras_enrichment dataset
        aminoacids = list(DEMO_DATASETS['df_sumo1'].index)

        # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_sumo1'].columns[0]

        # Full sequence
        sequence_sumo1 = 'MSDQEAKPSTEDLGDKKEGEYIKLKVIGQDSSEIHFKVKMTTHLKKLKESYCQRQGVPMN' + 'SLRFLFEGQRIADNHTPKELGMEEEDVIEVYQEQTGGHSTV'
        # Define secondary structure
        secondary_sumo1 = [['L0'] * (20), ['β1'] * (28 - 20), ['L1'] * 3, ['β2'] * (39 - 31),
                           ['L2'] * 4, ['α1'] * (55 - 43),
                           ['L3'] * (6), ['β3'] * (65 - 61), ['L4'] * (75 - 65), ['α2'] * (80 - 75),
                           ['L5'] * (85 - 80), ['β4'] * (92 - 85), ['L6'] * (101 - 92)]

        return Screen(
            DEMO_DATASETS['df_sumo1'], sequence_sumo1, aminoacids, start_position, 1,
            secondary_sumo1
        )

    def _generate_mapk1_obj(self) -> Screen:
        """
        Create object for the mapk1 dataset.
        """
        # Order of amino acid substitutions in the hras_enrichment dataset
        aminoacids = list(DEMO_DATASETS['df_mapk1'].index)

        # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_mapk1'].columns[0]

        # Full sequence
        sequence_mapk1_x = 'MAAAAAAGAGPEMVRGQVFDVGPRYTNLSYIGEGAYGMVCSAYDNVNKVRVAIK' + 'KISPFEHQTYCQRTLREIKILLRFRHENIIGINDIIRAPTIEQMKDVYIVQDLMETDLYKLLKTQ' + 'HLSNDHICYFLYQILRGLKYIHSANVLHRDLKPSNLLLNTTCDLKICDFGLARVADPDHDHTGFL' + 'TEYVATRWYRAPEIMLNSKGYTKSIDIWSVGCILAEMLSNRPIFPGKHYLDQLNHILGILGSPSQ' + 'EDLNCIINLKARNYLLSLPHKNKVPWNRLFPNADSKALDLLDKMLTFNPHKRIEVEQALAHPYLE' + 'QYYDPSDEPIAEAPFKFDMELDDLPKEKLKELIFEETARFQPGYRS'

        # Create objects
        return Screen(DEMO_DATASETS['df_mapk1'], sequence_mapk1_x, aminoacids, start_position, 0)

    def _generate_ube2i_obj(self) -> Screen:
        """
        Create object for the ube2i dataset.
        """

        # Order of amino acid substitutions in the hras_enrichment dataset
        aminoacids = list(DEMO_DATASETS['df_ube2i'].index)

        # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_ube2i'].columns[0]

        # Full sequence
        sequence_ube2i_x = 'MSGIALSRLAQERKAWRKDHPFGFVAVPTKNPDGTMNLMNWECAIPGKKGTP' + 'WEGGLFKLRMLFKDDYPSSPPKCKFEPPLFHPNVYPSGTVCLSILEEDKDWRPAITIKQ' + 'ILLGIQELLNEPNIQDPAQAEAYTIYCQNRVEYEKRVRAQAKKFAPS'

        # Define secondary structure
        secondary_ube2i = [['α1'] * (20 - 1), ['L1'] * (24 - 20), ['β1'] * (30 - 24), ['L2'] * 5,
                           ['β2'] * (46 - 35), ['L3'] * (56 - 46), ['β3'] * (63 - 56),
                           ['L4'] * (73 - 63), ['β4'] * (77 - 73), ['L5'] * (93 - 77),
                           ['α2'] * (98 - 93), ['L6'] * (107 - 98), ['α3'] * (122 - 107),
                           ['L7'] * (129 - 122), ['α4'] * (155 - 129), ['L8'] * (160 - 155)]

        # Create objects
        return Screen(
            DEMO_DATASETS['df_ube2i'], sequence_ube2i_x, aminoacids, start_position, 1,
            secondary_ube2i
        )

    def _generate_tat_obj(self) -> Screen:
        """
        Create object for the tat dataset.
        """

        # Order of amino acid substitutions in the hras_enrichment dataset
        aminoacids = list(DEMO_DATASETS['df_tat'].index)

        # First residue of the hras_enrichment dataset. Because 1-Met was not mutated, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_tat'].columns[0]

        # Full sequence
        sequence_tat = 'MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKALGISYGRKKRRQRRRAHQ' + 'NSQTHQASLSKQPTSQPRGDPTGPKE'

        # Define secondary structure
        secondary_tat = [['L1'] * (8), ['α1'] * (13 - 8), ['L2'] * (28 - 14), ['α2'] * (41 - 28),
                         ['L3'] * (90 - 41)]

        return Screen(
            DEMO_DATASETS['df_tat'], sequence_tat, aminoacids, start_position, 0, secondary_tat
        )

    def _generate_rev_obj(self) -> Screen:
        """
        Create object for the rev dataset.
        """

        # Order of amino acid substitutions in the hras_enrichment dataset
        aminoacids = list(DEMO_DATASETS['df_rev'].index)

        # First residue of the hras_enrichment dataset. Because 1-Met was not mureved, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_rev'].columns[0]

        # Full sequence
        sequence_rev = 'MAGRSGDSDEDLLKAVRLIKFLYQSNPPPNPEGTRQARRNRRRRWRERQRQIHSISERIL' + 'STYLGRSAEPVPLQLPPLERLTLDCNEDCGTSGTQGVGSPQILVESPTILESGAKE'

        # Define secondary structure
        secondary_rev = [['L1'] * (8), ['α1'] * (25 - 8), ['L2'] * (33 - 25), ['α2'] * (68 - 33),
                         ['L3'] * (116 - 41)]

        return Screen(
            DEMO_DATASETS['df_rev'], sequence_rev, aminoacids, start_position, 0, secondary_rev
        )

    def _generate_asynuclein_obj(self) -> Screen:
        """
        Create object for the synuclein dataset.
        """

        # Order of amino acid substitutions in the hras_enrichment dataset
        aminoacids = list(DEMO_DATASETS['df_asynuclein'].index)

        # First residue of the hras_enrichment dataset. Because 1-Met was not mureved, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_asynuclein'].columns[0]

        # Full sequence
        sequence_asynuclein = 'MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTK' + 'EQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDP' + 'DNEAYEMPSEEGYQDYEPEA'

        # Define secondary structure
        secondary_asynuclein = [['L1'] * (1), ['α1'] * (37 - 1), ['L2'] * (44 - 37),
                                ['α2'] * (92 - 44), ['L3'] * (140 - 92)]

        return Screen(
            DEMO_DATASETS['df_asynuclein'], sequence_asynuclein, aminoacids, start_position, 0,
            secondary_asynuclein
        )

    def _generate_aph_obj(self) -> Screen:
        """
        Create object for the aph dataset.
        """

        aminoacids = list(DEMO_DATASETS['df_aph'].index)

        # First residue of the hras_enrichment dataset. Because 1-Met was not mureved, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_aph'].columns[0]

        # Full sequence
        sequence_aph = 'MIEQDGLHAGSPAAWVERLFGYDWAQQTIGCSDAAVFRLSAQGRPVLFVKTDLSGALNELQ' + 'DEAARLSWLATTGVPCAAVLDVVTEAGRDWLLLGEVPGQDLLSSHLAPAEKVSIMADAMRR' + 'LHTLDPATCPFDHQAKHRIERARTRMEAGLVDQDDLDEEHQGLAPAELFARLKARMPDGED' + 'LVVTHGDACLPNIMVENGRFSGFIDCGRLGVADRYQDIALATRDIAEELGGEWADRFLVLY' + 'GIAAPDSQRIAFYRLLDEFF'

        # Define secondary structure
        secondary_aph = [['L1'] * (11), ['α1'] * (16 - 11), ['L2'] * (22 - 16), ['β1'] * (26 - 22),
                         ['L3'] * (34 - 26), ['β2'] * (40 - 34), ['L4'] * (46 - 40), ['β3'] *
                         (52 - 46), ['L5'] * (58 - 52), ['α2'] * (72 - 58), ['L6'] * (79 - 72),
                         ['β4'] * (85 - 79), ['L7'] * (89 - 85), ['β5'] * (95 - 89),
                         ['L8'] * (99 - 95), ['β6'] * (101 - 99), ['L9'] * (107 - 101),
                         ['α3'] * (131 - 107), ['L10'] * (135 - 131), ['α4'] * (150 - 135),
                         ['L11'] * (158 - 150), ['α5'] * (163 - 158), ['L12'] * (165 - 163),
                         ['α6'] * (177 - 165), ['L13'] * (183 - 177), ['β7'] * (187 - 183),
                         ['L14'] * (191 - 187), ['α7'] * (194 - 191), ['L15'] * (1),
                         ['β8'] * (199 - 195), ['L16'] * (201 - 199), ['β9'] * (206 - 201),
                         ['L17'] * (212 - 206), ['β10'] * (216 - 212), ['α8'] * (245 - 216),
                         ['L18'] * (4), ['α9'] * (264 - 249)]

        return Screen(
            np.log10(DEMO_DATASETS['df_aph']), sequence_aph, aminoacids, start_position, 0,
            secondary_aph
        )

    def _generate_b11L5F_obj(self) -> Screen:
        """
        Create object for the aph dataset.
        """

        # Order of amino acid substitutions in the hras_enrichment dataset
        aminoacids = list(DEMO_DATASETS['df_b11L5F'].index)

        # Sequence
        sequence_b11L5F = 'CRAASLLPGTWQVTMTNEDGQTSQGQMHFQPRSPYTLDVKAQGTISDGRPI' + 'SGKGKVTCKTPDTMDVDITYPSLGNMKVQGQVTLDSPTQFKFDVTTSDGSKVTGTLQRQE'

        # First residue of the hras_enrichment dataset. Because 1-Met was not mureved, the dataset starts at residue 2
        start_position = DEMO_DATASETS['df_b11L5F'].columns[0]

        return Screen(DEMO_DATASETS['df_b11L5F'], sequence_b11L5F, aminoacids, start_position, 0)
