"""
This module contains the set up instructions.
"""
import sys
import os
import setuptools
from mutagenesis_visualization import (__author__, __author_email__, __version__,
                                        __title__, __license__, __description__)

if sys.version_info[:2] < (3, 8):
    raise RuntimeError("Python version >= 3.8 required.")

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
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics'],
    project_urls={
        'Documentation': ('https://mutagenesis-visualization.readthedocs.io/'),
        'Publication': 'https://doi.org/10.1101/2021.10.08.463725',
        'Source': 'https://github.com/fhidalgor/mutagenesis_visualization',},
    install_requires=read_list('requirements.txt'),
    dependency_links=['https://github.com/cxhernandez/ipymol/tarball/master'],
    zip_safe=True,
    include_package_data=True,
    python_requires='>=3.8',
    )
