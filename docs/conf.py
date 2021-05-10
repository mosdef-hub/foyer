# -*- coding: utf-8 -*-
#
# Configuration file for the Sphinx documentation builder.
#
# This file does only contain a selection of the most common options. For a
# full list see the documentation:
# http://www.sphinx-doc.org/en/stable/config

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#


import os
import sys
import pathlib

sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('sphinxext'))

base_path = pathlib.Path(__file__).parent
os.system('python {} --name'.format((base_path / '../setup.py').resolve()))


import foyer


# -- Project information -----------------------------------------------------

project = u'foyer'
copyright = u'2014-2020, Vanderbilt University'

# The short X.Y version
version = '0.8.1'
# The full version, including alpha/beta/rc tags
release = version


# -- General configuration ---------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#
# needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'nbsphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.autosummary'
]

# nbsphinx specific configuration
nbsphinx_execute = 'never'
nbsphinx_input_prompt = '%s'
nbsphinx_prompt_width = 1.1
html_scaled_image_link = False
nbsphinx_allow_errors = False



# Prolog: Displayed on top of the notebook

nbsphinx_prolog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}
.. raw:: html
  <div class="admonition note">
      <p>This page was generated from
        <a class="reference external" href="https://github.com/mosdef-hub/foyer/blob/{{ env.config.release|e }}/{{ docname|e }}">{{ docname|e }}</a>.
        <br>Interactive online version:
        <a href="https://mybinder.org/v2/gh/mosdef-hub/foyer_tutorials/master?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a>.
      </p>
    </div>
.. raw:: latex
    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from
    \sphinxcode{\sphinxupquote{\strut {{ docname | escape_latex }}}} \dotfill}}
"""

# Displayed on bottom on notebook
nbsphinx_epilog = r"""
{% set docname = 'docs/' + env.doc2path(env.docname, base=None) %}
.. raw:: latex

\nbsphinxstopnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{\dotfill\ \sphinxcode{\sphinxupquote{\strut
    {{ docname | escape_latex }}}} ends here.}}
"""


autosummary_generate = True
autodoc_default_flags = ['members', 'inherited-members']
numpydoc_class_members_toctree = False

napoleon_numpy_docstring = True
napoleon_use_admonition_for_notes = True
napoleon_use_param = False


_python_doc_base = 'http://docs.python.org/3.7'


intersphinx_mapping = {
    _python_doc_base: None,
    'http://docs.scipy.org/doc/numpy': None,
    'http://docs.scipy.org/doc/scipy/reference': None,
    'http://scikit-learn.org/stable': None
}

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
#language = None

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path .
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'sphinx'


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'
import sphinx_rtd_theme
html_theme = 'sphinx_rtd_theme'
hhtml_theme_path= [
    sphinx_rtd_theme.get_html_theme_path()
]

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
#
# html_theme_options = {}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# Custom sidebar templates, must be a dictionary that maps document names
# to template names.
#
# The default sidebars (for documents that don't match any pattern) are
# defined by theme itself.  Builtin themes are using these templates by
# default: ``['localtoc.html', 'relations.html', 'sourcelink.html',
# 'searchbox.html']``.
#
html_sidebars = {
'**': ['globaltoc.html', 'sourcelink.html', 'searchbox.html', 'relations.html', 'sourcelink.html'],
}

# -- Options for HTMLHelp output ---------------------------------------------

# If true, the index is split into individual pages for each
html_split_index = False

# Output file base name for HTML help builder.
htmlhelp_basename = 'foyerdoc'


# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    #
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    #
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #
    # 'preamble': '',

    # Latex figure (float) alignment
    #
    # 'figure_align': 'htbp',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title,
#  author, documentclass [howto, manual, or own class]).
latex_documents = [
    ('index', 'foyer.tex', u'foyer Documentation',
     u'Author', 'manual'),
]


# -- Options for manual page output ------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'foyer', u'foyer Documentation',
     [u'Author'], 1)
]


# -- Options for Texinfo output ----------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
    ('index', 'foyer', u'foyer Documentation',
     u'Author', 'foyer', 'One line description of project.',
     'Miscellaneous'),
]

# -- Options for Epub output ----------------------------------------------

# Bibliographic Dublin Core info.
epub_title = u'foyer'
epub_author = u'Author'
epub_publisher = u'Author'
epub_copyright = u'2020, Author'

# The basename for the epub file. It defaults to the project name.
#epub_basename = u'foyer'

# The HTML theme for the epub output. Since the default themes are not optimized
# for small screen space, using the same theme for HTML and epub output is
# usually not wise. This defaults to 'epub', a theme designed to save visual
# space.
#epub_theme = 'epub'

# The language of the text. It defaults to the language option
# or en if the language is not set.
#epub_language = ''

# The scheme of the identifier. Typical schemes are ISBN or URL.
#epub_scheme = ''

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#epub_identifier = ''

# A unique identification for the text.
#epub_uid = ''

# A tuple containing the cover image and cover page html template filenames.
#epub_cover = ()

# A sequence of (type, uri, title) tuples for the guide element of content.opf.
#epub_guide = ()

# HTML files that should be inserted before the pages created by sphinx.
# The format is a list of tuples containing the path and title.
#epub_pre_files = []

# HTML files shat should be inserted after the pages created by sphinx.
# The format is a list of tuples containing the path and title.
#epub_post_files = []

# A list of files that should not be packed into the epub file.
epub_exclude_files = ['search.html']

# The depth of the table of contents in toc.ncx.
#epub_tocdepth = 3

# Allow duplicate toc entries.
#epub_tocdup = True

# Choose between 'default' and 'includehidden'.
#epub_tocscope = 'default'

# Fix unsupported image types using the PIL.
#epub_fix_images = False

# Scale large images.
#epub_max_image_width = 0

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#epub_show_urls = 'inline'

# If false, no index is generated.
#epub_use_index = True


# -- Extension configuration -------------------------------------------------
