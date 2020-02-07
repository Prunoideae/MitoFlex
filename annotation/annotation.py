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
import json

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call, concat_command
    from utility.profiler import profiling
    # TODO remove annotation when testing
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


def annotate(basedir=None, prefix=None, ident=30, fastafile=None,
             genetic_code=9, clade=None, taxa=None, thread_number=8,
             wildcard_profile=False, trna_overlapping=40):
    logger.log(2, 'Entering annotation module.')
    if wildcard_profile:
        logger.log(
            3, 'Wildcard protein profile is used, results may not be accurate.')

    # Once we can confirm the sequences are from the clade we want to,
    # then we don't need to use overall database.
    if wildcard_profile:
        # Copying the code because I'm lazy here
        logger.log(2, 'Updating the general protein database.')
        lc = 0
        with open(path.join(profile_dir_tbn, 'Animal.fa'), 'w') as fout:
            for protein_fas in os.listdir(profile_dir_tbn):
                if protein_fas.endswith('.fa') and protein_fas != 'Animal.fa':
                    with open(path.join(profile_dir_tbn, protein_fas)) as fin:
                        for line in fin:
                            fout.write(line)
                            lc += 1
        logger.log(1, f'Generation finished with {lc} writes.')

    tbn_profile = path.join(
        profile_dir_tbn, f'{clade if not wildcard_profile else "Animal"}.fa')
    blast_file = tk.tblastn(dbfile=tbn_profile, infile=fastafile, genetic_code=genetic_code,
                            basedir=basedir, prefix=prefix)

    blast_frame, _ = tk.blast_to_csv(blast_file, ident=ident, score=25)
    blast_frame = tk.wash_blast_results(blast_frame)

    wise_frame, _, _ = tk.genewise(
        basedir=basedir, prefix=prefix, wises=blast_frame,
        infile=fastafile, dbfile=tbn_profile, cutoff=0.5)

    # Output wise frame
    wise_csv = path.join(basedir, f'{prefix}.genewise.result.csv')
    wise_frame.to_csv(wise_csv)
    logger.log(1, f'Genewise results generated at {wise_csv}')

    trna_out_dir = path.join(basedir, 'trna')
    os.makedirs(trna_out_dir, exist_ok=True)
    query_dict, missing_trna = tk.trna_search(
        fastafile, profile_dir_trna, trna_out_dir, prefix, genetic_code, 0.01, overlap_cutoff=trna_overlapping)

    logger.log(2, f'tRNAs found : {list(query_dict.keys())}')
    if missing_trna:
        logger.log(3, f'Missing tRNAs : {missing_trna}')

    rrna_out_dir = path.join(basedir, 'rrna')
    os.makedirs(rrna_out_dir, exist_ok=True)
    result_12, result_16 = tk.rrna_search(
        fastafile, profile_dir_rrna, rrna_out_dir, prefix, 0.01)

    if not result_12:
        logger.log(3, '12s rRNA is not found!')

    if not result_16:
        logger.log(3, '16s rRNA is not found!')

    locs_file = path.join(basedir, 'locs.json')
    annotation_json = {}

    sequence_data = {x.id: x for x in SeqIO.parse(fastafile, 'fasta')}

    annotated_fa = path.join(basedir, f'{prefix}.annotated.cds.fa')
    annotated_frag = []
    for _, row in wise_frame.iterrows():
        cds = str(row).split('_')[3]
        if cds in annotation_json:
            count = sum(x.startswith(cds) for x in annotation_json.keys())
            cds = f'{cds}{"_" if count > 0 else ""}{count}'
            if count > 0:
                logger.log(
                    3, f'Duplicated gene {cds} detected at {start} - {end}!')
        start, end = (min(int(row.wise_min_start), int(row.wise_max_end)),
                      max(int(row.wise_min_start), int(row.wise_max_end)))
        frag = sequence_data[str(row.sseq)][start-1:end]
        frag.description = f'gene={cds} start={start} end={end}'
        annotated_frag.append(frag)
        annotation_json[cds] = (start, end, 0, str(row.sseq))

    SeqIO.write(annotated_frag, annotated_fa, 'fasta')

    annotated_rnas = path.join(basedir, f'{prefix}.annotated.rna.fa')
    annotated_frag.clear()
    for key, value in query_dict.items():
        start, end = (min(value.seqfrom, value.seqto),
                      max(value.seqfrom, value.seqto))
        frag = sequence_data[value.sequence][start-1:end]
        frag.description = f'gene=trn{key} start={start} end={end}'
        annotated_frag.append(frag)
        annotation_json[f'trn{key}'] = (
            start, end, 1, value.sequence)

    if result_12:
        start, end = (min(result_12.seqfrom, result_12.seqto),
                      max(result_12.seqfrom, result_12.seqto))
        logger.log(
            2, f'12s rRNA found from {start} to {end}')
        frag = sequence_data[result_12.sequence][start-1:end]
        frag.description = f'gene=rrnS start={start} end={end}'
        annotated_frag.append(frag)
        annotation_json['rrnS'] = (
            start, end, 2, result_12.sequence)

    if result_16:
        start, end = (min(result_16.seqfrom, result_16.seqto),
                      max(result_16.seqfrom, result_16.seqto))
        logger.log(
            2, f'16s rRNA found from {start} to {end}')
        frag = sequence_data[result_16.sequence][start-1:end]
        frag.description = f'gene=rrnL start={start} end={end}'
        annotated_frag.append(frag)
        annotation_json['rrnL'] = (
            start, end, 2, result_16.sequence)

    SeqIO.write(annotated_frag, annotated_rnas, 'fasta')
    with open(locs_file, 'w') as f:
        json.dump(annotation_json, f)

    return locs_file, annotated_fa, annotated_rnas
