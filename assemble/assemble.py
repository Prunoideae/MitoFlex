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
    from utility import logger
    from configurations import assemble as a_conf  # Prevent naming confliction
    from assemble.assemble_wrapper import MEGAHIT, EmptyGraph  # pylint: disable=import-error, no-name-in-module
    from assemble.scaffold_wrapper import SOAP, scaf2mega
except ImportError as err:
    sys.exit(
        f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)


def assemble(fastq1=None, fastq2=None, base_dir=None, work_prefix=None,
             kmer_list=None, depth_list=None, disable_local=False,
             prune_level=2, prune_depth=2, keep_temp=False,
             threads=8, min_multi=3.0, insert_size=125, no_scaf=False):

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
        'min_depth': min_multi,
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
        if not disable_local:
            megahit.local(c, n)
        megahit.iterate(c, n)

    megahit.finalize(megahit.kmax)

    if not no_scaf and fastq2 is not None:
        soap = SOAP(fastq1, fastq2, megahit.final_contig,
                    libread.max_len, insert_size, base_dir, threads, work_prefix, megahit.kmax)
        logger.log(2, "Building lib.")
        soap.lib()
        logger.log(2, "Calling SOAP-Wrapper.")
        return soap.scaf()
    elif fastq2 is None:
        logger.log(2, "Scaffolding skipped due to SE reads.")
    else:
        logger.log(2, "Scaffolding skipped due to disabled.")

    return megahit.final_contig
