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
import os
import sys
import datetime
sys.path.append(os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../../'))
import xgrads


# -- Project information -----------------------------------------------------

project = 'xgrads'
copyright = f'{datetime.datetime.today().year}, MiniUFO'
author = 'MiniUFO'

# The short X.Y version
version = xgrads.__version__
# The full version, including alpha/beta/rc tags
release = xgrads.__version__


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'sphinx.ext.autodoc',     # api auto-gen
    'sphinx.ext.doctest',
    'sphinx.ext.todo',
    'sphinx.ext.mathjax',     # math
    # 'myst_parser',            # Markdown
    'recommonmark',           # Markdown
]

# The master toctree document.
# master_doc = 'index'

# The suffix(es) of source filenames, either a string or list.
source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = [
    'conf.py', 'sphinxext', '_build', '_templates', '_themes',
    '*.ipynb', '**.ipynb_checkpoints' '.DS_Store', 'trash', 'tmp',
]


# -- Options for HTML output -------------------------------------------------

# Logo
html_logo = os.path.join('_static', 'xgradsLogo.png')

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': True,
    'display_version': False,
    'collapse_navigation': True,
    'navigation_depth': 4,
    'prev_next_buttons_location': 'bottom',  # top and bottom
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large. Static folder is for CSS and image files. Use ImageMagick to
# convert png to ico on command line with 'convert image.png image.ico'
html_favicon = os.path.join('_static', 'xgradsIcon.ico')


