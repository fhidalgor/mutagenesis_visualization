from os import system
from mutagenesis_visualization import __version__

system("poetry run jupyter nbconvert --ClearOutputPreprocessor.enabled=True --clear-output --inplace mutagenesis_visualization/tutorial/doc*.ipynb")
system("poetry run jupyter nbconvert --output-dir='docs/' --to rst mutagenesis_visualization/tutorial/doc*.ipynb")

# edit conf.py version
with open("docs/conf.py", "r", encoding="utf8") as f:
    lines = f.readlines()
with open("docs/conf.py", "w", encoding="utf8") as f:
    for line in lines:
        if line.startswith("version = "):
            f.write('version: "{}" '.format(__version__))
        else:
            f.write(line)

# edit image loading from doc0_intro.rst
with open("docs/doc0_intro.rst", "r", encoding="utf8") as f:
    lines = f.readlines()
with open("docs/doc0_intro.rst", "w", encoding="utf8") as f:
    for line in lines:
        if line.strip("\n") == ".. |image0| image:: ../../docs/_static/workflow_v3.png":
            f.write(".. |image0| image:: _static/workflow_v3.png\n")
        else:
            f.write(line)

# edit version from CITATION.cff
with open("CITATION.cff", "r", encoding="utf8") as f:
    lines = f.readlines()
with open("CITATION.cff", "w", encoding="utf8") as f:
    for line in lines:
        if line.startswith("version: "):
            f.write('version: "{}" '.format(__version__))
        else:
            f.write(line)

system("poetry export --without-hashes > requirements.txt")
system("git add .")
system("git commit -m 'minor change to parameter'")
system("git push")
