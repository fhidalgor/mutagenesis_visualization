poetry run yapf -ir mutagenesis_visualization/
poetry run mypy --config-file mypy.ini mutagenesis_visualization/
poetry run pylint -d duplicate-code mutagenesis_visualization/
