from distutils.core import setup, Extension
try:
   from distutils.command.build_py import build_py_2to3 as build_py
except ImportError:
   from distutils.command.build_py import build_py

import sys, os
try:
    import argparse
except ImportError:
    sys.exit('This program requires Python 2.7.2 or higher, or Python 3.2 or higher. Please upgrade your version (%s) of Python and try again.' % (sys.version.split()[0]))


setup(name = "pheno_man",
    version = '0.1.0',
    description = "Sample manipulation program for ESP phenotype data",
    author = 'Biao Li and Gao Wang',
    author_email = 'biaol@bcm.edu',
    maintainer = 'Biao Li',
    maintainer_email = 'lb4@rice.edu',
    py_modules = [
        'pheno_man.__init__',
        'pheno_man.pheno_man',
        'pheno_man.str4r'
    ],
    scripts = ['phenoman'],
    cmdclass = {'build_py': build_py },
    package_dir = {'pheno_man': 'pheno_man'},
    packages = ['pheno_man'],
    package_data = {'pheno_man': ['EA.txt', 'AA.txt', 'EAdups.txt', 'AAdups.txt', 'ALL.txt']}
)
