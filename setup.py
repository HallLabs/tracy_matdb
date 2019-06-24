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
      version='1.5.0',
      description='Database generator for ML materials space.',
      long_description= "" if not path.isfile("README.md") else read_md('README.md'),
      author='HALL LABS',
      author_email='wmorgan@tracy.com',
      url='https://github.com/NewVistas/tracy_matdb',
      license='GNU GPLv3',
      setup_requires=['pytest-runner',],
      tests_require=['pytest',
                     'numpy',
                     'phonopy',
                     'termcolor',
                     'aflux'],
      install_requires=[
          "numpy",
          "argparse",
          "ase",
          "pyparsing",
          "termcolor",
          "six",
          "numpy",
          "requests",
          "beautifulsoup4",
          "tqdm",
          "html5lib",
          "mpld3",
          "phenum",
          "h5py",
          "seekpath",
          "PyYAML"
      ],
      packages=['matdb', 'matdb.database', 'matdb.fitting','matdb.calculators','matdb.plotting'],
      scripts=['matdb/scripts/matdb_build.py',
               'matdb/scripts/matdb_plot.py',
               'matdb/scripts/matdb_train.py',
               'matdb/scripts/matdb_find.py',
               'matdb/scripts/matdb_move.py',
               'matdb/scripts/matdb_supercell.py',
               'matdb/scripts/matdb_convert.py',
               'matdb/scripts/matdb_mtp_to_relax.py',
               'matdb/scripts/matdb_fix.py',
               'matdb/scripts/matdb_watch.py',
               'support/matdb_vasp.py',
               'support/matdb_sbatch.py',
               'support/matdb_module.py',
               'support/matdb_mlp.py',
               'support/matdb_getkpoints.py'],
      package_dir={
          'matdb': 'matdb'
      },
      package_data={
          'matdb': [
              'matdb/templates/*',
              'matdb/templates/uniqueBinaries/*',
              'matdb/templates/uniqueTernaries/*',
              'matdb/templates/uniqueUnaries/*',
      ]},
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
