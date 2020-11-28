# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
import sys
sys.path.append('../')
#sys.path.insert(0, os.path.abspath('../mutagenesis_visualization/'))
#sys.path.insert(0, os.path.abspath('../../mutagenesis_visualization/'))

# -- Project information -----------------------------------------------------
project = 'Mutagenesis Visualization'
copyright = '2020, Frank Hidalgo'
author = 'Frank Hidalgo'

# The full version, including alpha/beta/rc tags
version = '0.0.2'


# -- General configuration ---------------------------------------------------
# If your documentation needs a minimal Sphinx version, state it here.
needs_sphinx = '1.7.3'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon', 'numpydoc',
'sphinx.ext.autosectionlabel', 'sphinx.ext.intersphinx', 'IPython.sphinxext.ipython_console_highlighting',
'IPython.sphinxext.ipython_directive']

# mock import modules
autodoc_mock_imports = ['Bio', 'ipymol', 'math', 'logomaker', 'collections', 'adjustText', 'seaborn',
'scipy', 'sklearn', 'pandas', 'Import_notebook','copy', 'itertools', 'freesasa', 'plotly', 'statsmodels']

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
#
# source_suffix = ['.rst', '.md']
source_suffix = '.rst'

# The master toctree document.
master_doc = 'index'

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#
# This is also used if you do content translation via gettext catalogs.
# Usually you set "language" from the command line for these cases.
language = None


# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.

html_theme_options = {
    'canonical_url': 'https://mutagenesis-visualization.readthedocs.io',
    'logo_only': True,
    'display_version': True,
    'prev_next_buttons_location': 'none', #'bottom',
    'style_external_links': False,
    #'vcs_pageview_mode': '',
    'style_nav_header_background': 'white',
    # Toc options
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 3,
    'includehidden': True,
    'titles_only': False
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

html_logo = '_static/logo.png'

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    ('index', 'main.tex', 'Mutagenesis Visualization Documentation',
     'Frank Hidalgo', 'manual')
 ]
     
latex_engine = 'pdflatex'

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
    'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
    'pointsize': '12pt',}

# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [('index', 'mutagenesis-visualization', 'Mutagenesis Visualization Documentation',
     ['Frank Hidalgo'], 1)]
     