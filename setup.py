import setuptools

def readme():
    with open("README.md", "r") as fh:
        return fh.read()

setuptools.setup(
    name="mutagenesis_visualization", # Replace with your own username
    version="0.0.1",
    author="Frank Hidalgo",
    author_email="fhidalgoruiz@berkeley.edu",
    description="Package for analysis and visualization of saturation mutagenesis data",
    long_description=readme(),
    url="https://github.com/fhidalgor/mutagenesis_visualization",
    keywords=['saturation','mutagenesis','deepsequencing'],
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3.0",
        "Operating System :: OS Independent",
    ],
    install_requires=['__future__', 'numpy', 'seaborn', 'pandas', 'copy', 'scipy', 'matplotlib',
        'sklearn', 'adjustText', 'itertools', 'ipymol', 'pymol', 'collections', 'Bio', 'collections',
        'logomaker'],
    zip_safe=False,
    include_package_data=True,
    python_requires='>=3.6',