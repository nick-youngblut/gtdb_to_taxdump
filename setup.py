#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import glob
import numpy
from bin import __version__

    
# dependencies
install_reqs = [
   'networkx'
]

## install main application
desc = 'GTDB database utility scripts'
setup(
    name = 'gtdb_to_taxdump',
    version = '0.1.7',
    description = desc,
    long_description = desc + '\n See README for more information.',
    author = 'Nick Youngblut',
    author_email = 'nyoungb2@gmail.com',
    install_requires = install_reqs,
    include_dirs = [numpy.get_include()],
    packages = find_packages(),
    package_dir={'gtdb2td' : 'gtdb2td'},
    license = "MIT license",
    url = 'https://github.com/nick-youngblut/gtdb_to_taxdump',
    scripts = ['bin/gtdb_to_diamond.py',
               'bin/gtdb_to_taxdump.py',
               'bin/lineage2taxid.py',
               'bin/ncbi-gtdb_map.py']
)




