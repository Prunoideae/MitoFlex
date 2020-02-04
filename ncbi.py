#!/usr/bin/env python3

"""
ncbi.py
========

Copyright (c) 2019-2020 Li Junyu <2018301050@szu.edu.cn>.

This file is part of MitoFlex.

MitoFlex is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoFlex is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoFlex.  If not, see <http://www.gnu.org/licenses/>.

"""

import os
from os.path import getmtime
from datetime import datetime

import sys
if sys.version_info[0] < 3:
    sys.exit('Python 3 must be installed in current environment! Please check if your environment setup (like conda environment) is deactivated or wrong!')

try:
    from ete3 import NCBITaxa
except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoFlex installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

if os.path.isfile('taxdump.tar.gz'):
    os.rename('taxdump.tar.gz', 'old.taxdump.tar.gz')
try:
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    if os.path.isfile('old.taxdump.tar.gz'):
        os.remove('old.taxdump.tar.gz')
except Exception as idd:
    print("Errors occured when fetching data from NCBI database, falling back to the last fetched database.")
    ncbi = NCBITaxa(taxdump_file=os.path.abspath('old.taxdump.tar.gz'))
    if os.path.isfile('old.taxdump.tar.gz'):
        os.rename('old.taxdump.tar.gz', 'taxdump.tar.gz')
