"""
annotation.py
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
import sys
from os import path
import subprocess
import multiprocessing

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call, concat_command
    from utility.profiler import profiling
    # TODO: remove annotation when testing
    from annotation import annotation_tookit as tk  # pylint: disable=import-error, no-name-in-module
    from Bio import SeqIO
    from utility import logger
    from utility.bio import infernal
    from utility.bio import wuss
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

mitoflex_dir = path.abspath(path.join(path.dirname(__file__), '..'))
profile_dir = path.join(mitoflex_dir, 'profile')
profile_dir_hmm = path.join(profile_dir, 'CDS_HMM')
profile_dir_tbn = path.join(profile_dir, 'MT_database')
profile_dir_rrna = path.join(profile_dir, 'rRNA_CM')
profile_dir_trna = path.join(profile_dir, 'tRNA_CM')


def concat_java(*args, **kwargs):
    return concat_command(*args, **kwargs).replace('--', '-')


@profiling
def annotate(basedir=None, prefix=None, ident=30, fastafile=None,
               genetic_code=9, clade=None, taxa=None, thread_number=8):

    # Once we can confirm the sequences are from the clade we want to,
    # then we don't need to use overall database.
    tbn_profile = path.join(profile_dir_tbn, f'{clade}_CDS_protein.fa')
    blast_file = tk.tblastn(dbfile=tbn_profile, infile=fastafile, genetic_code=genetic_code,
                            basedir=basedir, prefix=prefix)

    blast_frame, _ = tk.blast_to_csv(blast_file, ident=ident, score=25)
    blast_frame = tk.wash_blast_results(blast_frame)

    wise_frame, queries, database = tk.genewise(
        basedir=basedir, prefix=prefix, wises=blast_frame,
        infile=fastafile, dbfile=tbn_profile, cutoff=0.5)

    # Output logging file
    final_file = path.join(basedir, f'{prefix}.genewise.predicted.cds.fa')
    if not tk.collect_result(final_file, wise_frame, queries, database):
        return None, None

    reloc = tk.reloc_genes(fasta_file=fastafile, wises=wise_frame)
    wise_csv = path.join(basedir, f'{prefix}.genewise.result.csv')
    reloc.to_csv(wise_csv)

    trna_out_dir = path.join(basedir, 'trna')
    os.makedirs(trna_out_dir, exist_ok=True)
    query_dict, missing_trna = tk.trna_search(
        fastafile, profile_dir_trna, trna_out_dir, prefix, genetic_code, 0.01)

    if missing_trna:
        logger.log(3, f'Missing tRNAs for AAs : {missing_trna}')

    rrna_out_dir = path.join(basedir, 'rrna')
    os.makedirs(rrna_out_dir, exist_ok=True)
    result_12, result_16 = tk.rrna_search(
        fastafile, profile_dir_rrna, rrna_out_dir, prefix, 0.01)

    if not result_12:
        logger.log(3, '12s rRNA is not found!')

    if not result_16:
        logger.log(3, '16s rRNA is not found!')

    return wise_csv, query_dict, result_12, result_16
