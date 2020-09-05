import setuptools

def readme():
    with open("README.md", "r") as fh:
        return fh.read()

setuptools.setup(
    name="mutagenesis_visualization",
    version="0.0.5",
    author="Frank Hidalgo",
    author_email="fhidalgoruiz@berkeley.edu",
    description="A package for analysis and visualization of saturation mutagenesis data",
    long_description=readme(),
    url="https://github.com/fhidalgor/mutagenesis_visualization",
    keywords=['saturation','mutagenesis','deepsequencing'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3.0",
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Education',
        'Natural Language :: English',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ],
    keywords=('ras', 'mutagenesis', 'saturation', 'saturation mutagenesis', 'deepsequencing'),
    project_urls={
        'Documentation': ('https://htmlpreview.github.io/?https://github.com/'
                          'drewhart/geonomics/blob/master/doc/built/doc.html'),
        # 'Methods Paper': 'PAPER URL HERE!',
        'Source': 'https://github.com/fhidalgor/mutagenesis_visualization',
    install_requires=['numpy>=1.18.5', 'seaborn>=0.10.1', 'pandas>=1.0.5', 'copy', 'scipy>=1.5.0', 'matplotlib>=3.2.2',
        'sklearn>=0.23.1', 'adjustText>=0.7.3', 'itertools>=8.4.0', 'ipymol>=0.5', 'collections>=1.2.1', 'Bio>=1.77']
        'logomaker>=0.8', 'Shannon>=1.0.0'],
    zip_safe=True,
    include_package_data=True,
    python_requires='>=3.6',