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
    from utility.helper import shell_call
    from utility import logger
    from configurations import assemble as a_conf  # Prevent naming confliction
except Exception:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")


def assemble(fastq1=None, fastq2=None, base_dir=None, work_prefix=None,
             uselist=False, kmin=21, kmax=141, kstep=12, klist=None,
             disable_local=False,
             prune_level=2, prune_depth=2, keep_temp=False,
             threads=8, addtional_kmers=[]):

    logger.log(2, 'Start assembling mitochondrial sequences.')

    if(uselist):
        kmin = kmax = kstep = None
        logger.log(1, f'Using kmer list : {klist}')
    else:
        klist = None
        logger.log(1, f'Using step parameters : min={kmin}, max={kmax}')

    logger.log(
        1, f'Using arguments : no_local={disable_local} ,p_lv = {prune_level}, p_dep = {prune_depth}')

    tmp_dir = path.join(base_dir, 'temp')
    try:
        os.makedirs(tmp_dir, exist_ok=True)
    except Exception:
        tmp_dir = None

    options = {
        'k_min': kmin,
        'k_max': kmax,
        'k_step': kstep,
        'k_list': klist,
        'no_mercy': a_conf.no_mercy,
        'prune_level': prune_level,
        'prune_depth': prune_depth,
        'keep_tmp_files': keep_temp,
        'tmp_dir': tmp_dir,
        'out_dir': path.join(base_dir, 'result'),
        'out_prefix': work_prefix,
        'no_hw_accel': a_conf.disable_acc,
        'num_cpu_threads': threads,
        'no_local': disable_local,
        'min_count': a_conf.min_multi,
        'kmin_1pass': a_conf.one_pass
    }

    logger.log(0, f'Calling megahit with : {options}')

    # Mutable options
    if fastq1 and fastq2:
        options['_1'] = fastq1
        options['_2'] = fastq2
    elif fastq1:
        options['r'] = fastq1
    elif fastq2:
        options['r'] = fastq2

    shell_call('megahit', **options)

    contigs_file = os.path.join(
        base_dir, 'result', work_prefix + '.contigs.fa')
    if not path.isfile(contigs_file):
        raise RuntimeError(
            'Megahit finished with no results! This could be an error caused memory overflow or something will halt the process of Megahit!')
    if logger.get_level() <= 1:
        from Bio import SeqIO
        contigs = [x for x in SeqIO.parse(contigs_file, 'fasta')]
        logger.log(
            1, f'Output contigs file size : {path.getsize(contigs_file)}')
        logger.log(1, f'Contig number : {len(contigs)}')

    inter_dir = path.join(base_dir, "result", "intermediate_contigs")

    if addtional_kmers:
        logger.log(
            1, f"Merging intermediate contigs {addtional_kmers} with finals.")

        shell_call(
            "cat",
            *[path.join(inter_dir, f'k{x}.contigs.fa')
              for x in addtional_kmers],
            '>>',
            contigs_file
        )

    if not keep_temp:
        logger.log(1, f'Cleaning intermidiate contig files.')
        shell_call("rm -r", inter_dir)

    return contigs_file


# This is currently NOT a working function.
# The effect of scaffolding using SOAPdenovo-fusion is still
# under investigation.
def scaffolding(fastq1=None, fastq2=None):
    pass
