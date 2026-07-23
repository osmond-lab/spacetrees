import os
import sys

sys.path.insert(0, os.path.abspath('..'))

project = 'spacetrees'
copyright = '2024, Osmond Lab'
author = 'Matthew Osmond'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'myst_parser',
]

autodoc_mock_imports = ['tskit', 'tsconvert']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
