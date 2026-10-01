import re
from pathlib import Path

from setuptools import setup, find_packages

# ...read out of the package rather than retyped here. Importing it would pull
# in numpy and IPython before they are installed, so this reads the line.
_init = Path(__file__).parent.joinpath('py2Dmol', '__init__.py').read_text()
VERSION = re.search(r'^__version__ = "([^"]+)"', _init, re.M).group(1)

setup(
    name='protein-inspector',
    version='1.13.0',
    # AUTHORSHIP OF THIS PACKAGE, not of the viewer it is built on. These
    # fields said 'sokrypton' / so3@mit.edu / the py2Dmol URL, which credited
    # generously but falsely: it implied Sergey Ovchinnikov wrote and endorsed
    # this derivative, and it would have sent its bug reports to him. He is
    # credited where credit is due -- LICENSE, NOTICE, and every exported
    # bundle -- and the upstream URL is below as a project URL rather than as
    # this package's home.
    author='Tadas Kluonis',
    description='Manifest-driven structural inspection bundles, built on py2Dmol.',
    long_description=(
        'Single-file interactive structure inspection bundles: annotation '
        'layers, conformation morphing and partner toggles over the py2Dmol '
        'viewer by Sergey Ovchinnikov, which is inlined into every export and '
        'credited there. Both parts are BEER-WARE (Revision 42). See NOTICE.'
    ),
    long_description_content_type='text/markdown',
    url='https://github.com/tadaskluonis/protein_inspector',
    project_urls={'Upstream viewer (py2Dmol)': 'https://github.com/sokrypton/py2Dmol'},
    # NOTICE requires both files to travel with the package, and without this
    # a built wheel shipped neither.
    license_files=['LICENSE', 'NOTICE'],
    packages=find_packages(),
    include_package_data=True,
    # EVERY RESOURCE viewer.py OPENS. It reads these by name through
    # importlib.resources, so a file missing here is not a degraded viewer - it
    # is a FileNotFoundError on the first show(), in the wheel only, where no
    # test in this repo runs. viewer-cartoon-gpu.min.js was missing and
    # viewer.py:1329 opens it unconditionally.
    #
    # There is no MANIFEST.in and no pyproject.toml, so include_package_data
    # above contributes nothing and this list is the whole of it.
    # tests/packaging.py builds a wheel and imports it, which is the only way
    # this list can be checked at all.
    # PER-DIRECTORY GLOBS, not eighteen literal paths. The resources moved into
    # core/, parts/, cartoon/, panels/ and align/, and a list that long is a list
    # that goes stale - which it did once already, shipping a wheel without the
    # GPU renderer. tests/packaging.py checks these against what viewer.py reads.
    # WHAT SHIPS IS THE BUNDLES, not the two dozen source files they are built
    # from. tools/bundle.py concatenates and minifies each target; the panels
    # stay loose because viewer.py adds them only when the config asks.
    # tests/packaging.py checks this against what viewer.py actually reads.
    package_data={
        'py2Dmol': [
            'resources/viewer.html',
            # ...the NOTEBOOK bundle only. The glob that was here also shipped
            # py2Dmol.embed.min.js, a web artefact in
            # every pip install, which viewer.py never opens. tests/packaging.py
            # fails if it ever opens one that is not listed here.
            # ONE, since the three narrow notebook builds became one complete
            # bundle - see tools/bundle.py. Naming files that no longer exist
            # is not an error setuptools reports; it just ships nothing for
            # them, and the omission only shows up in a release environment
            # where the setuptools-scm plugin is not there to cover for it.
            'resources/bundles/py2Dmol.notebook.min.js',
        ],
    },
    license='BEER-WARE',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    # LOWER BOUNDS, NOT PINS. `==` in a library's install_requires is a
    # declaration that it will not co-exist with anything, and pip resolves it
    # against every other package in the user's environment. The pins here were
    # whatever happened to be installed the day the file was written; numpy
    # 2.4.6 alone already implies Python >= 3.11, which `>=3.6` denied.
    python_requires='>=3.10',
    install_requires=[
        'numpy>=1.24',
        'ipython>=8.0',
        'biopython>=1.81',
        # NOT OPTIONAL, whatever the docs used to say. py2Dmol/viewer.py does
        # `import gemmi` at module scope (line ~213), so it is needed to import
        # this package at all, not merely to read a .pdb. Listing it as an
        # extra made `pip install protein-inspector` produce an installation
        # that raises ModuleNotFoundError on the first import -- reproducible
        # only from a built wheel in a clean environment, which is exactly
        # where no test in this repo used to look.
        'gemmi>=0.6',
    ],
    extras_require={
        # The MCP server, for agent hosts that speak the protocol.
        'mcp': ['mcp>=1.2'],
        # Taking the picture without a human pressing Save. A browser draws
        # it, so this pulls a browser driver in; the bundle itself never
        # needs one.
        'capture': ['playwright>=1.40'],
    },
    entry_points={
        'console_scripts': [
            'protein-inspector = protein_inspector.__main__:cli',
            'protein-inspector-mcp = protein_inspector.mcp_server:main',
        ],
    },
)
