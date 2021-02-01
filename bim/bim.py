"""
bim.py
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
import sys
from os import path
from typing import Tuple, Union

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call

except ImportError as err:
    sys.exit(
        f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)


def bim_assemble(threads: int, fasta_file: str, prefix: str, fastq1: str, fastq2: Union[str, None] = None) -> Tuple[str, Union[str, None]]:
    '''
    A special assemble function that chains tons of program together to assemble from bait and reads.
    '''
