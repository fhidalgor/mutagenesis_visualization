"""
This module contains the set up instructions.
"""
import sys
import os
import setuptools
from mutagenesis_visualization import (__author__, __author_email__, __version__,
                                        __title__, __license__, __description__)

if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python version >= 3.6 required.")

def readme():
    """returns readme file."""
    with open("README.md", "r", enconding="utf8") as fh:
        return fh.read()

def read(filename: str):
    """read file"""
    base_dir = os.path.dirname(__file__)
    filename = os.path.join(base_dir, filename)
    with open(filename, 'r', enconding="utf8") as fi:
        return fi.read()

def read_list(filename: str):
    """reads requirements.txt"""
    rows = read(filename).split("\n")
    rows = [x.strip() for x in rows if x.strip()]
    return list(rows)


setuptools.setup(
    name=__title__,
    version=__version__,
    author=__author__,
    author_email=__author_email__,
    description= __description__,
    long_description=readme(),
    url="https://github.com/fhidalgor/mutagenesis_visualization",
    keywords=['saturation','mutagenesis','deepsequencing', 'ras', 'site-saturation', 'nextgenseq', 'NGS'],
    packages=setuptools.find_packages(),
    #py_modules=['main.mutagenesis_visualization'],
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    project_urls={
        'Documentation': ('https://mutagenesis-visualization.readthedocs.io/'),
        # 'Methods Paper': 'PAPER URL HERE!',
        'Source': 'https://github.com/fhidalgor/mutagenesis_visualization',},
#    install_requires=['numpy>=1.19.5', 'seaborn>=0.10.0', 'pandas>=1.2.0', 'scipy>=1.5.0', 'matplotlib>=3.3',
#        'scikit-learn>=0.24.2', 'adjustText>=0.7.3', 'biopython>=1.79','freesasa>=2.1.0', 'plotly>=5.1.0'
#        'statsmodels>=0.12.2', 'xlsxwriter>=1.4.4', 'xlrd>=2.0.1'],
    install_requires=read_list('requirements.txt'),
    dependency_links=['https://github.com/cxhernandez/ipymol/tarball/master'],
    zip_safe=True,
    include_package_data=True,
    python_requires='>=3.6',
    )
