"""
In this module I will combine all the data paths so it is easier to
keep track and maintain.
"""
import pkg_resources

PATH_VARIABLE = "mutagenesis_visualization"

# PDBs
PDB_5P21 = pkg_resources.resource_filename(PATH_VARIABLE, "data/5p21.pdb")
PDB_1ERM = pkg_resources.resource_filename(PATH_VARIABLE, "data/1erm.pdb")
PDB_1A5R = pkg_resources.resource_filename(PATH_VARIABLE, "data/1a5r.pdb")
PDB_1ND4 = pkg_resources.resource_filename(PATH_VARIABLE, "data/1nd4.pdb")

# XLSX
HRAS_RBD_COUNTS = pkg_resources.resource_filename(PATH_VARIABLE, "data/hrasRBD_counts.xlsx")
HRAS_GAPGEF_COUNTS = pkg_resources.resource_filename(PATH_VARIABLE, "data/hrasGAPGEF_counts.xlsx")

# CSV
HRAS_RBD_COUNTS_CSV = pkg_resources.resource_filename(PATH_VARIABLE, "data/HRas166_RBD.csv")
HRAS_GAPGEF_COUNTS_CSV = pkg_resources.resource_filename(PATH_VARIABLE, "data/HRas166_GAPGEF.csv")

# PKL
DF_BLA_RAW_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_bla_raw.pkl")
DF_SUMO1_RAW_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_sumo1_raw.pkl")
DF_MAPK1_RAW_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_mapk1_raw.pkl")
DF_UBE2I_RAW_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_ube2i_raw.pkl")
DF_TAT_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_tat.pkl")
DF_REV_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_rev.pkl")
DF_ASYNUCLEIN_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_asynuclein.pkl")
DF_APH_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_aph.pkl")
DF_B11L5F_RAW_PKL = pkg_resources.resource_filename(PATH_VARIABLE, "data/df_b11l5f_raw.pkl")

# FASTA and FASTQ
DEMO_FASTA = pkg_resources.resource_filename(PATH_VARIABLE, "data/Ras_family_trimmed.fasta")
HRAS_FASTQ = pkg_resources.resource_filename(PATH_VARIABLE, "data/hras.trimmed.fastq")
TEST_SHORT_FASTQ = pkg_resources.resource_filename(PATH_VARIABLE, "data/for_tests/short.fastq")

# LOGO
PATH_LOGO = pkg_resources.resource_filename(PATH_VARIABLE, 'data/logo.xlsx')
