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
import json
import warnings

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import concat_command
    from annotation import annotation_tookit as tk  # pylint: disable=import-error, no-name-in-module
    from Bio import SeqIO
    from Bio import BiopythonWarning
    from utility import logger
    from utility.bio import seq
    from misc.check_circular import check_circular
    import configurations
except ImportError as identifier:
    sys.exit(
        f"Unable to import helper module {identifier.name}, is the installation of MitoFlex valid?")

mitoflex_dir = path.abspath(path.join(path.dirname(__file__), '..'))
profile_dir = path.join(mitoflex_dir, 'profile')
profile_dir_hmm = path.join(profile_dir, 'CDS_HMM')
profile_dir_tbn = path.join(profile_dir, 'MT_database')
profile_dir_rrna = path.join(profile_dir, 'rRNA_CM')
profile_dir_trna = path.join(profile_dir, 'tRNA_CM')


def concat_java(*args, **kwargs):
    return concat_command(*args, **kwargs).replace('--', '-')


def annotate(basedir=None, prefix=None, ident=30, fastafile=None,
             genetic_code=9, clade=None, thread_number=8,
             wildcard_profile=False, trna_overlapping=40, hmmer_search=True, score=5, e_value=0.005):
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

    blast_file = tk.tblastn_multi(dbfile=tbn_profile, infile=fastafile, genetic_code=genetic_code,
                                  basedir=basedir, prefix=prefix, threads=thread_number)
    blast_frame, _ = tk.blast_to_csv(blast_file, ident=ident, score=25)
    washed_frame = tk.wash_blast_results(blast_frame)

    if configurations.annotation.redirection:
        logger.log(2, 'Checking genome directions.')
        if tk.redirect_genome(fasta_file=fastafile, blast_frame=blast_frame):
            # Not fixing frame directly because I'm lazy, since later I'll need to fix tRNA and rRNA frames if so...
            logger.log(2, "Genome is reversed, launching a second annotation to fix gene locations.")
            blast_file = tk.tblastn_multi(dbfile=tbn_profile, infile=fastafile, genetic_code=genetic_code,
                                          basedir=basedir, prefix=prefix, threads=thread_number)
            blast_frame, _ = tk.blast_to_csv(blast_file, ident=ident, score=25)
            washed_frame = tk.wash_blast_results(blast_frame)

    wise_frame, _, _ = tk.genewise(
        basedir=basedir, prefix=prefix, wises=washed_frame,
        infile=fastafile, dbfile=tbn_profile, cutoff=0.5)

    # Add an extra washing here, since Pandas will have some
    # strange behaviour processing data this large...
    # Also prevent some mutation of ['plus']
    wise_frame = tk.wash_blast_results(wise_frame, mut_plus=False)

    taxa_data = {}
    for _, row in wise_frame.iterrows():
        splited = str(row.qseq).split('_')
        PCG = splited[3]
        taxa_name = ' '.join(splited[4:6])
        taxa_score = float(row.score)

        if PCG not in taxa_data or taxa_data[PCG][1] < taxa_score:
            taxa_data[PCG] = (taxa_name, taxa_score)

    score_data = {}
    for k, (taxa_name, taxa_score) in taxa_data.items():
        score_data[taxa_name] = score_data[taxa_name] + \
            1 if taxa_name in score_data else 0

    most_possible = list(score_data.keys())[0]
    for taxa_name, taxa_score in score_data.items():
        if score_data[most_possible] < taxa_score:
            most_possible = taxa_name

    logger.log(2, f'Determined most possible species : {most_possible}')

    if configurations.annotation.reloc_genes:
        logger.log(2, 'Relocating genes.')
        wise_frame = tk.reloc_genes(fasta_file=fastafile,
                                    wises=wise_frame, code=genetic_code)

    cds_indexes = {}
    cds_found = []
    with open(path.join(profile_dir_hmm, 'required_cds.json')) as f:
        cds_indexes = json.load(f)[clade]

    for _, row in wise_frame.iterrows():
        cds = str(row).split('_')[3]
        cds_found.append(cds)

    hmmer_frame = None

    cds_notfound = [x for x in cds_indexes if x not in cds_found]
    logger.log(2, f'PCGs found in annotation : {cds_found}')
    if cds_notfound and not hmmer_search:
        logger.log(3, f'Expected PCG {cds_notfound} not found!')
    elif cds_notfound and hmmer_search:
        logger.log(
            3, f'Expected PCG {cds_notfound} not found, turning to nhmmer search.'
        )
        hmmer_frame = tk.nhmmer_search(fasta_file=fastafile, thread_number=thread_number,
                                       nhmmer_profile=profile_dir_hmm + f'/{clade}.hmm', prefix=prefix, basedir=basedir)
        hmmer_frame = hmmer_frame[~hmmer_frame['query'].isin(cds_found)]
        hmmer_frame = hmmer_frame[hmmer_frame['e'] < e_value]
        hmmer_frame = hmmer_frame[hmmer_frame['score'] > score]
        logger.log(2, 'Recovered pcgs : \n' + str(hmmer_frame))

    trna_out_dir = path.join(basedir, 'trna')
    os.makedirs(trna_out_dir, exist_ok=True)

    # Disable some annoying warning
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
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
    start = end = -1
    for _, row in wise_frame.iterrows():
        cds = str(row).split('_')[3]
        if cds in annotation_json:
            count = sum(x.startswith(cds) for x in annotation_json.keys())
            cds = f'{cds}{"_" if count > 0 else ""}{count}'
        start, end = (min(int(row.wise_min_start), int(row.wise_max_end)),
                      max(int(row.wise_min_start), int(row.wise_max_end)))
        frag = sequence_data[str(row.sseq)][start - 1:end]
        frag.description = f'gene={cds} start={start} end={end} from={row.sseq} strand={"+" if row.plus else "-"}'
        annotated_frag.append(frag)
        annotation_json[cds] = (start, end, 0, str(
            row.sseq), "+" if row.plus else "-")

    if hmmer_frame is not None:
        for _, row in hmmer_frame.iterrows():
            start, end = (min(int(row.envfrom), int(row.envto)),
                          max(int(row.envfrom), int(row.envto)))
            frag = sequence_data[str(row.target)][start - 1:end]
            frag.description = f'gene={str(row.query)} start={start} end={end} from={row.target} strand={row.strand}'
            annotated_frag.append(frag)
            annotation_json[str(row.query)] = (
                start, end, 0, str(row.target), str(row.strand))

    SeqIO.write(annotated_frag, annotated_fa, 'fasta')

    annotated_rnas = path.join(basedir, f'{prefix}.annotated.rna.fa')
    annotated_frag.clear()
    for key, value in query_dict.items():
        start, end = (min(value.seqfrom, value.seqto),
                      max(value.seqfrom, value.seqto))
        frag = sequence_data[value.sequence][start - 1:end]
        frag.description = f'gene=trn{key} start={start} end={end}'
        annotated_frag.append(frag)
        annotation_json[f'trn{key}'] = (
            start, end, 1, value.sequence, '+' if value.plus else '-')

    if result_12:
        start, end = (min(result_12.seqfrom, result_12.seqto),
                      max(result_12.seqfrom, result_12.seqto))
        logger.log(
            2, f'12s rRNA found from {start} to {end}')
        frag = sequence_data[result_12.sequence][start - 1:end]
        frag.description = f'gene=rrnS start={start} end={end}'
        annotated_frag.append(frag)
        annotation_json['rrnS'] = (
            start, end, 2, result_12.sequence, '+' if result_12.plus else '-')

    if result_16:
        start, end = (min(result_16.seqfrom, result_16.seqto),
                      max(result_16.seqfrom, result_16.seqto))
        logger.log(
            2, f'16s rRNA found from {start} to {end}')
        frag = sequence_data[result_16.sequence][start - 1:end]
        frag.description = f'gene=rrnL start={start} end={end}'
        annotated_frag.append(frag)
        annotation_json['rrnL'] = (
            start, end, 2, result_16.sequence, '+' if result_16.plus else '-')

    SeqIO.write(annotated_frag, annotated_rnas, 'fasta')
    with open(locs_file, 'w') as f:
        json.dump(annotation_json, f, indent=4, separators=(',', ': '))

    return locs_file, annotated_fa, annotated_rnas


def fix_circular(fa_file: str):
    genome = [x for x in SeqIO.parse(fa_file, 'fasta')]
    circular = False
    # Only one sequence
    if len(genome) == 1:
        traits = seq.decompile(genome[0].description, sep=' ')
        # The sequence is circular
        if 'flag' in traits and traits['flag'] == '3':
            results = check_circular(final_fasta=fa_file)
            # The overlapped region is determined
            if results:
                result = results[0]
                if result[0] != -1:
                    # Trim the genome and returns
                    circular = True
                    logger.log(
                        2, f'An overlapped region was found starting at {result[0][0]} with length {result[0][1]}. Trimming it.')
                    end = genome[0].seq.rfind(result[1])
                    genome[0] = genome[0][result[0][0]:end]
                    with open(fa_file, 'w') as f:
                        SeqIO.write(genome, f, 'fasta')
    return circular
