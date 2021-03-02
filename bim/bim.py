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
    from utility.helper import direct_call
    from utility import logger
    from utility.helper import timed
except ImportError as err:
    sys.exit(
        f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)


@timed(enabled=True)
def bwa_map(threads: int, fasta_file: str, basedir: str, prefix: str,
            fastq1: str, fastq2: str, quality: int = 30) -> Tuple[str, str, str]:
    index = path.join(basedir, prefix)
    direct_call(f'bwa index -p {index} {fasta_file}')
    fq1, fq2 = path.join(basedir, prefix + '.1.fq'), path.join(basedir, prefix + '.2.fq') if fastq2 is not None else None
    bam = path.join(basedir, prefix + ".bam")
    logger.log(2, "Mapping and extracting reads from bwa mem.")
    direct_call(
        f'\
        bwa mem -t {threads} {index} {fastq1} {fastq2 if fastq2 is not None else ""} |\
        samtools view -bS -q {quality} -h - |\
        tee {bam}|\
        samtools fastq -1 {fq1} {f"-2 {fq2}" if fq2 is not None else ""} -')

    return bam, fq1, fq2


def cal_insert(bam: str, threads: int, basedir: str, prefix: str) -> int:
    stat_file = path.join(basedir, prefix + ".stats")
    stats = direct_call(
        f'\
        samtools stats {bam}|\
        tee {stat_file}|\
        grep ^IS|\
        cut -f 2-'
    )
    logger.log(2, str(stats))
    return 150
