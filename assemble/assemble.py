"""
assemble.py
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

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call
    from utility.profiler import profiling
    from utility import logger
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")


def assemble(fastq1=None, fastq2=None, base_dir=None, work_prefix=None,
             uselist=False, kmin=21, kmax=141, kstep=12, klist=None,
             no_mercy=False, disable_acc=False, disable_local=False,
             prune_level=2, prune_depth=2, keep_temp=False,
             threads=8):

    logger.log(2, 'Start assembling mitochondrial sequences.')

    if(uselist):
        kmin = kmax = kstep = None
        logger.log(1, f'Using step list : {klist}')
    else:
        klist = None
        logger.log(1, f'Using step parameters : min={kmin}, max={kmax}')

    logger.log(
        1, f'Using arguments : mercy={not no_mercy}, no_acc={disable_acc}, no_local={disable_local} ,p_lv = {prune_level}, p_dep = {prune_depth}')

    tmp_dir = path.join(base_dir, 'temp')
    try:
        os.makedirs(tmp_dir, exist_ok=True)
    except:
        tmp_dir = None

    shell_call('megahit', _1=fastq1, _2=fastq2,
               k_min=kmin, k_max=kmax, k_step=kstep, k_list=klist,
               no_mercy=no_mercy, prune_level=prune_level, prune_depth=prune_depth,
               keep_tmp_files=keep_temp, tmp_dir=tmp_dir,
               out_dir=path.join(base_dir, 'result'), out_prefix=work_prefix,
               no_hw_accel=disable_acc, num_cpu_threads=threads, no_local=disable_local)

    contigs_file = os.path.join(
        base_dir, 'result', work_prefix + '.contigs.fa')

    if logger.get_level() <= 1:
        from Bio import SeqIO
        contigs = [x for x in SeqIO.parse(contigs_file, 'fasta')]
        logger.log(
            1, f'Output contigs file size : {path.getsize(contigs_file)}')
        logger.log(1, f'Contig number : {len(contigs)}')

    return contigs_file
