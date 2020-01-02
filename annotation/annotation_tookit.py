"""
annotation_toolkit.py
=========

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
from os import path
import sys
try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import concat_command, direct_call, shell_call
    from utility.profiler import profiling
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)

# Truncates all the -- to - to suit blast's parsing style


def truncated_call(*args, **kwargs):
    direct_call(concat_command(*args, **kwargs).replace('--', '-'))


def filter_blast():
    pass


def tblastn(dbfile=None, infile=None, genetic_code=9, basedir=None,
            prefix=None):
    # Since 'in' is the preserve keyword
    truncated_call('makeblastdb', '-in', dbfile, dbtype='nucl')

    # Go tblastn
    out_blast = path.join(path.abspath(basedir), f'{prefix}.blast')
    truncated_call('tblastn', useconv=False, evalue='1e-5', outfmt=6,
                   seg='no', db_gencode=genetic_code, db=infile,
                   query=dbfile)
    out_blast_filtered = filter_blast()
    solar(out_blast_filtered)


def solar(blastfile=None):
    solar_path = path.join(bin_dir, 'solar', 'solar.pl')
    truncated_call('perl', solar_path, blastfile, '>', f'{blastfile}.solar')
