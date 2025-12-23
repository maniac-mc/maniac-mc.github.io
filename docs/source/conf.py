# Configuration file for the Sphinx documentation builder.

project = 'MANIAC-MC'
copyright = '2025, Simon Gravelle'
author = 'Simon Gravelle'
release = 'v0.4.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.mathjax',
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
html_theme = 'furo'
html_static_path = ['_static']
