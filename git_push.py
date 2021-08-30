from os import system

system("poetry run jupyter nbconvert --ClearOutputPreprocessor.enabled=True --clear-output --inplace mutagenesis_visualization/tutorial/doc*.ipynb")
system("poetry run jupyter nbconvert --output-dir='docs/' --to rst mutagenesis_visualization/tutorial/doc*.ipynb")

# edit image loading from doc0_intro.rst
with open("docs/doc0_intro.rst", "r") as f:
    lines = f.readlines()
with open("docs/doc0_intro.rst", "w") as f:
    for line in lines:
        if line.strip("\n") == ".. |image0| image:: ../../docs/_static/workflow_v3.png":
            f.write(".. |image0| image:: _static/workflow_v3.png\n")
        else:
            f.write(line)

system("poetry export --without-hashes > requirements.txt")
system("git add .")
system("git commit -m 'documentation tweak'")
system("git push")
