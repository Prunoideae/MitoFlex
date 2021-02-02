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
import subprocess
import sys
from os import path
from typing import Tuple, Union

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call
    from bim.piped_wrapper import MEGAHIT

except ImportError as err:
    sys.exit(
        f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)


def bim_assemble(threads: int, fasta_file: str, basedir: str, prefix: str,
                 fastq1: str, fastq2: Union[str, None] = None, disable_local=False,
                 prune_level=2, prune_depth=2, keep_temp=False, insert_size=125,
                 no_scaf=False
                 ) -> Tuple[str, Union[str, None]]:
    '''
    A special assemble function that chains tons of program together to assemble from bait and reads.
    '''

    index = path.join(basedir, prefix)
    direct_call(f'bwa index -p {index} {fasta_file}')
    fq1, fq2 = path.join(basedir, prefix + '.1.fq'), path.join(basedir, prefix + '.2.fq') if fastq2 is not None else None
    os.mkfifo(fq1)
    if fq2 is not None:
        os.mkfifo(fq2)

    feed = subprocess.Popen(
        f'\
        bwa mem -t {threads} {fastq1} {fastq2 if fastq2 is not None else ""} |\
        samtools view -S -q 30 -h - |\
        {path.join(bin_dir,"sp")} sam2fq --fastq1 {fq1} {f"--fastq2 {fq2}" if fq2 is not None else ""}')

    options = {
        'prune_level': prune_level,
        'prune_depth': prune_depth,
        'keep_temp': keep_temp,
        'basedir': basedir,
        'prefix': prefix,
        'threads': threads,
        'no_local': disable_local,
        'fq1': fq1,
        'fq2': fq2,
    }

    megahit = MEGAHIT(**options)
    megahit.build_lib()
    if feed.wait() != 0:
        pass
