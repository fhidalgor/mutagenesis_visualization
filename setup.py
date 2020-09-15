import setuptools
import sys

if sys.version_info[:2] < (3, 6):
    raise RuntimeError("Python version >= 3.6 required.")

def readme():
    with open("README.md", "r") as fh:
        return fh.read()

setuptools.setup(
    name="mutagenesis_visualization",
    version="0.0.1.4",
    author="Frank Hidalgo",
    author_email="fhidalgoruiz@berkeley.edu",
    description="A package for analysis and visualization of saturation mutagenesis data",
    long_description=readme(),
    url="https://github.com/fhidalgor/mutagenesis_visualization",
    keywords=['saturation','mutagenesis','deepsequencing', 'ras'],
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
    install_requires=['numpy', 'seaborn', 'pandas', 'scipy', 'matplotlib>=3.0.0',
        'sklearn', 'adjustText', 'Bio', 'logomaker','freesasa'],
    dependency_links=['https://github.com/cxhernandez/ipymol.git'],
    zip_safe=True,
    include_package_data=True,
    python_requires='>=3.6',
    )