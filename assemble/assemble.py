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
    from assemble.assemble_wrapper import MEGAHIT  # pylint: disable=import-error, no-name-in-module
except Exception:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)


def assemble(fastq1=None, fastq2=None, base_dir=None, work_prefix=None,
             kmer_list=None, depth_list=None, disable_local=False,
             prune_level=2, prune_depth=2, keep_temp=False,
             threads=8, min_multi=3.0):

    # TODO : Modify the megahit's default workflow to make it optimized.

    logger.log(2, 'Start assembling mitochondrial sequences.')

    logger.log(1, f'Using kmer list : {kmer_list}')

    options = {
        'prune_level': prune_level,
        'prune_depth': prune_depth,
        'keep_temp': keep_temp,
        'basedir': base_dir,
        'prefix': work_prefix,
        'threads': threads,
        'no_local': disable_local,
        'fq1': fastq1,
        'fq2': fastq2,
    }

    logger.log(2, f'Initializing megahit wrapper.')
    megahit = MEGAHIT(**options)

    megahit.initialize()
    libread = megahit.build_lib()

    if libread.max_len + 20 < kmer_list[-1]:
        logger.log(
            3, f'Input max read length {libread.max_len} < max k-mer length {kmer_list[-1]}, resizing.')
        kmer_list = [*[x for x in kmer_list if x < libread.max_len + 20],
                     libread.max_len + 20]
        logger.log(3, f'K-mers after resized : {kmer_list}')

    megahit.kmax = kmer_list[-1]

    kmer_list = [0, *kmer_list, -1]
    kmer_list = [(kmer_list[i],
                  kmer_list[i + 1],
                  kmer_list[i + 2])
                 for i in range(len(kmer_list) - 2)]

    for i, (p, c, n) in enumerate(kmer_list):
        megahit.graph(p, c)
        megahit.assemble(c)
        megahit.filter(c, min_depth=depth_list[i],
                       min_length=0 if n != -1 else a_conf.min_length, max_length=a_conf.max_length)
        if n == -1:
            break
        megahit.local(c, n)
        megahit.iterate(c, n)

    megahit.finalize()

    return megahit.final_contig


# This is currently NOT a working function.
# The effect of scaffolding using SOAPdenovo-fusion is still
# under investigation.
def scaffolding(fastq1=None, fastq2=None):
    soap_fusion = path.join(bin_dir, 'SOAPdenovo-fusion')
    soap_127 = path.join(bin_dir, 'SOAPdenovo-127mer')
    pass
