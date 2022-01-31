"""
In this module the most important functions/classes will be loaded.
"""
from mutagenesis_visualization.main.classes.screen import Screen
from mutagenesis_visualization.main.classes.counts import Counts
from mutagenesis_visualization.main.classes.create_variants import CreateVariants
from mutagenesis_visualization.main.classes.generate_primers import GeneratePrimers
from mutagenesis_visualization.main.demo.demo_figures import run_demo
from mutagenesis_visualization.main.demo.demo_objects import DemoObjects
from mutagenesis_visualization.main.demo.demo_data import load_demo_datasets
from mutagenesis_visualization.main.process_data.calculate_enrichment import calculate_enrichment
from mutagenesis_visualization.main.process_data.count_fastq import count_fastq
from mutagenesis_visualization.main.process_data.count_reads import count_reads

__author__ = "Frank Hidalgo"
__version__ = "1.0.0"
__title__ = "mutagenesis_visualization"
__license__ = "GPLv3"
__author_email__ = "fhidalgoruiz@berkeley.edu"
__description__ = "A package for processing, analysis and visualization of site-saturation mutagenesis data."
