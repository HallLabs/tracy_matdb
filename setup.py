#!/usr/bin/env python
try:
    from setuptools import setup
    args = {}
except ImportError:
    from distutils.core import setup
    print("""\
*** WARNING: setuptools is not found.  Using distutils...
""")

from setuptools import setup
try:
    from pypandoc import convert
    read_md = lambda f: convert(f, 'rst')
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    read_md = lambda f: open(f, 'r').read()

from os import path
setup(name='matdb',
      version='1.2.0',
      description='Database generator for ML materials space.',
      long_description= "" if not path.isfile("README.md") else read_md('README.md'),
      author='Conrad W Rosenbrock',
      author_email='rosenbrockc@gmail.com',
      url='https://github.com/rosenbrockc/matdb',
      license='MIT',
      setup_requires=['pytest-runner',],
      tests_require=['pytest', 'numpy', 'phonopy'],
      install_requires=[
          "argparse",
          "ase",
          "pyparsing",
          "termcolor",
          "six",
          "numpy",
          "phonopy",
          "requests",
          "beautifulsoup4",
          "tqdm",
          "html5lib",
          "mpld3",
          "phenum",
          "h5py",
          "lazy_import",
          "seekpath"
      ],
      packages=['matdb', 'matdb.database', 'matdb.fitting','matdb.calculators'],
      scripts=['matdb/scripts/matdb_build.py',
               'matdb/scripts/matdb_plot.py',
               'matdb/scripts/matdb_train.py',
               'matdb/scripts/matdb_find.py',
               'matdb/scripts/tracy_sub.py',
               'matdb/scripts/matdb_move.py',
               'matdb/scripts/matdb_supercell.py',
               'support/matdb_vasp.py',
               'support/matdb_sbatch.py',
               'support/matdb_module.py',
               'support/matdb_mlp.py',
               'support/matdb_getkpoints.py'],
      package_data={'matdb': []},
      include_package_data=True,
      classifiers=[
          'Development Status :: 4 - Beta',
          'Intended Audience :: Science/Research',
          'Intended Audience :: Developers',
          'Natural Language :: English',
          'Operating System :: MacOS',
          'Operating System :: Microsoft :: Windows',
          'Programming Language :: Python',
          'Programming Language :: Python :: 2',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.4',
          'Topic :: Scientific/Engineering',
      ],
     )
