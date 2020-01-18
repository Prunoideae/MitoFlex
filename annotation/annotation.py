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
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

mitoflex_dir = path.abspath(path.join(path.dirname(__file__), '..'))
profile_dir = path.join(mitoflex_dir, 'profile')
profile_dir_hmm = path.join(profile_dir, 'CDS_HMM')
profile_dir_tbn = path.join(profile_dir, 'MT_database')
profile_dir_rna = path.join(profile_dir, 'rRNA_CM')


def concat_java(*args, **kwargs):
    return concat_command(*args, **kwargs).replace('--', '-')


@profiling
def annotation(basedir=None, prefix=None, ident=30, fastafile=None,
               genetic_code=9, clade=None, taxa=None, thread_number=8):
    # Once we can confirm the sequences are from the clade we want to,
    # then we don't need to use overall database.
    tbn_profile = path.join(profile_dir_tbn, f'{clade}_CDS_protein.fa')
    blast_file = tk.tblastn(dbfile=tbn_profile, infile=fastafile, genetic_code=genetic_code,
                            basedir=basedir, prefix=prefix, ident=ident)

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

    mitfi_bin_dir = path.join(mitoflex_dir, 'annotation', 'mitfi')
    mitfi = 'mitfi.jar'

    # Prepare single file
    split_scaf_dir = path.join(basedir, 'split_scaf')
    os.mkdir(split_scaf_dir)
    for seq in SeqIO.parse(fastafile, 'fasta'):
        SeqIO.write(seq, path.join(split_scaf_dir,
                                   f'{seq.id}.splited.fa'), 'fasta')

    # trna_file = path.join(basedir, f'{prefix}.trna')
    tasks = []
    results = []
    for f in [path.join(split_scaf_dir, x) for x in os.listdir(split_scaf_dir) if x.endswith('fa')]:
        tasks.append(concat_java('java', f'-Xmx2048m',
                                 jar=path.join(mitfi_bin_dir, mitfi),
                                 cores=1, code=genetic_code, evalue=0.001, onlycutoff=True,
                                 appending=[f]))
    pool = multiprocessing.Pool(thread_number)

    def work(command):
        return direct_call(command)

    res = pool.map_async(work, tasks, callback=results.append)
    res.wait()
    print(results)
    # TODO: Analyze results
