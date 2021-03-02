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
    from utility.helper import direct_call
    from assemble.assemble_wrapper import MEGAHIT, EmptyGraph
    from utility import logger
    from configurations import assemble as a_conf
    from assemble.scaffold_wrapper import SOAP
    from utility.helper import timed
except ImportError as err:
    sys.exit(
        f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)


@timed(True)
def bim_assemble(threads: int, fasta_file: str, basedir: str, prefix: str,
                 fastq1: str, fastq2: Union[str, None] = None, disable_local=False,
                 prune_level=2, prune_depth=2, keep_temp=False, insert_size=125,
                 no_scaf=False, kmer_list=None, depth_list=None
                 ) -> str:
    '''
    A special assemble function that chains tons of program together to assemble from bait and reads.
    '''

    index = path.join(basedir, prefix)
    direct_call(f'bwa index -p {index} {fasta_file}')
    fq1, fq2 = path.join(basedir, prefix + '.1.fq'), path.join(basedir, prefix + '.2.fq') if fastq2 is not None else None

    logger.log(2, "Mapping and extracting reads from bwa mem.")
    direct_call(
        f'\
        bwa mem -t {threads} {index} {fastq1} {fastq2 if fastq2 is not None else ""} |\
        samtools view -bS -q 30 -h - |\
        samtools fastq -1 {fq1} {f"-2 {fq2}" if fq2 is not None else ""} -')

    options = {
        'prune_level': prune_level,
        'prune_depth': prune_depth,
        'keep_temp': keep_temp,
        'basedir': basedir,
        'prefix': prefix,
        'threads': threads,
        'no_local': disable_local,
        'min_depth': 3.0,
        'fq1': fq1,
        'fq2': fq2,
    }

    logger.log(2, f'Initializing megahit wrapper.')
    megahit = MEGAHIT(**options)

    megahit.initialize()
    libread = megahit.build_lib()

    logger.log(1, f"Loaded {libread.read_count} reads from bwa mem.")

    if libread.max_len + 20 < kmer_list[-1]:
        logger.log(
            3, f'Input max read length {libread.max_len} < max k-mer length {kmer_list[-1]}, resizing.')
        kmer_list = [x for x in kmer_list if x < libread.max_len + 20]
        if libread.max_len % 2 != 0:
            kmer_list.append(libread.max_len + 20)
        logger.log(3, f'K-mers after resized : {kmer_list}')

    megahit.kmax = kmer_list[-1]
    megahit.kmin = kmer_list[0]

    kmer_list = [0, *kmer_list, -1]
    kmer_list = [(kmer_list[i],
                  kmer_list[i + 1],
                  kmer_list[i + 2])
                 for i in range(len(kmer_list) - 2)]

    for i, (p, c, n) in enumerate(kmer_list):
        try:
            megahit.graph(p, c)
        except EmptyGraph:
            logger.log(
                3, f'Iteration broke at kmer = {p}, since no valid contig in kmer = {c} is done!')
            # Return to last iteration kmer sets.
            megahit.kmax = kmer_list[i - 1][0]
            break

        contig_info, _ = megahit.assemble(c)
        contig_filtered, *_ = megahit.filter(c, min_depth=depth_list[i], force_filter=c == megahit.kmax,
                                             min_length=0 if n != -1 else a_conf.min_length, max_length=a_conf.max_length,
                                             deny_number=a_conf.filter_keep)
        logger.log(
            1, f'Contig for kmer = {c} : {contig_filtered}/{contig_info.count}')

        if n == -1:
            break
        megahit.local(c, n)
        megahit.iterate(c, n)

    megahit.finalize(megahit.kmax)

    if not no_scaf:
        soap = SOAP(fastq1, fastq2, megahit.final_contig,
                    libread.max_len, insert_size, basedir, threads, prefix, megahit.kmax)
        logger.log(2, "Building lib.")
        soap.lib()
        logger.log(2, "Calling SOAP-Wrapper.")
        return soap.scaf()
    else:
        logger.log(2, "Scaffolding skipped due to disabled.")

    os.remove(fq1)
    if fq2 is not None:
        os.remove(fq2)

    return megahit.final_contig
