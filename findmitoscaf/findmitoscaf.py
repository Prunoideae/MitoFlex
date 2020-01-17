"""
findmitoscaf.py
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
import json
import operator

import pandas
import numpy as np
from Bio import SeqIO
from ete3 import NCBITaxa
from os import path


try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call, maxs
    from utility.profiler import profiling
    from utility.seq import compile_seq, decompile
    from annotation import annotation_tookit as tk
except Exception:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

ncbi = NCBITaxa()
mitoflex_dir = path.abspath(path.join(path.dirname(__file__), '..'))
profile_dir = path.join(mitoflex_dir, 'profile')
profile_dir_hmm = path.join(profile_dir, 'CDS_HMM')
profile_dir_tbn = path.join(profile_dir, 'MT_database')
profile_dir_rna = path.join(profile_dir, 'rRNA_CM')

rank_list = ['kindom', 'phylum', 'class',
             'order', 'family', 'genus', 'species']


def get_rank(taxa_name=None):
    name_dict = ncbi.get_name_translator([taxa_name])

    rank_dict = {
        'kindom': 'NA',
        'phylum': 'NA',
        'class': 'NA',
        'order': 'NA',
        'family': 'NA',
        'genus': 'NA',
        'species': 'NA'
    }

    for taxid in ncbi.get_lineage(name_dict[taxa_name][0]):
        rank = ncbi.get_rank(taxid)[taxid]
        taxa = ncbi.get_taxid_translator([taxid])[taxid]
        if rank in rank_dict:
            rank_dict[rank] = taxa

    return [(tax_class, tax_id) for tax_class, tax_id in rank_dict.items()]


@profiling
def findmitoscaf(thread_number=8, clade=None, prefix=None,
                 basedir=None, gene_code=9, taxa=None,
                 contigs_file=None, relaxing=0, multi=10, cover_valve=1, min_multi=3.0):

    # Drop all the sequences where multi is too low to do further analysis
    filtered_fa = f'{prefix}.contigs.filtered.fa'
    filtered_contigs = []
    for seq in SeqIO.parse(contigs_file, 'fasta'):
        trait_string = seq.description.replace(seq.id + ' ', '', 1)
        trait = decompile(trait_string, sep=' ')
        if float(trait['multi']) >= min_multi:
            filtered_contigs.append(seq)
    SeqIO.write(filtered_contigs, filtered_fa, 'fasta')

    # Do nhmmer search and collect, filter results
    nhmmer_profile = path.join(profile_dir_hmm, f'{clade}_CDS.hmm')

    # do hmmer search
    hmm_frame = nhmmer_search(fasta_file=contigs_file, thread_number=thread_number,
                              nhmmer_profile=nhmmer_profile, prefix=prefix,
                              basedir=basedir)

    # filter by taxanomy
    if taxa is not None:
        # We use an overall protein dataset to determine what clades diffrent seqs belonged to.
        tbn_profile = path.join(profile_dir_tbn, f'Animal_CDS_protein.fa')
        hmm_frame = filter_taxanomy(
            taxa=taxa, fasta_file=contigs_file, hmm_frame=hmm_frame,
            basedir=basedir, prefix=prefix, dbfile=tbn_profile, gene_code=gene_code,
            relaxing=relaxing)

    contig_data = {x.id: x
                   for x in SeqIO.parse(contigs_file, 'fasta')
                   if x.id in hmm_frame.target}

    # filter by multi
    contig_data_high = []
    contig_data_low = []

    for contig in contig_data:
        if contig.description.startswith(contig.id + ' '):
            contig.description = contig.description.replace(
                contig.id + ' ', '', 1)
        traits = decompile(contig.description, sep=' ')
        if float(traits['multi']) >= multi:
            # Append traits to avoid parsing again
            contig_data_high.append((contig, traits))
        else:
            contig_data_low.append(contig)
            # Here we dispose all the low abundance contigs,
            # so only hmm_frame and contigs_file_high will be used.
            hmm_frame = hmm_frame[hmm_frame.target != contig.id]

    contigs_file_high = path.join(basedir, f'{prefix}.abundance.high.fa')
    contigs_file_low = path.join(basedir, f'{prefix}.abundance.low.fa')

    SeqIO.write([x[0] for x in contig_data_high], contigs_file_high, 'fasta')
    SeqIO.write(contig_data_low, contigs_file_low, 'fasta')

    # Here we pick out the last sequences by using the brute
    import brute

    cds_indexes = list(
        json.load(open(path.join(profile_dir, 'required_cds.json')))[clade])

    # First, collect the matching range of single sseq
    seq_dict = {}
    by_seqid = dict(tuple(hmm_frame.groupby['target']))
    for key, frame in by_seqid.items():
        found_cds = [(0, 0)] * len(cds_indexes)
        for _, row in frame.iterrows():
            query = str(row.query)
            query_start = int(row.hmmstart)
            query_end = int(row.hmmend)
            found_cds[cds_indexes.index(query)] = (query_start, query_end)
        seq_dict[found_cds] = key

    # Then, solve the problem by brute
    best = brute.solution(list(seq_dict.keys()))

    # Finally, convert returned solution to cds
    picked_ids = [seq_dict[seq_id] for seq_id in best]
    picked_seq = [seq for seq in contig_data if seq.id in picked_ids]
    missing_pcgs = [cds_indexes[zero] for zero in best if zero == (0, 0)]
    found_pcgs = [
        others for others in cds_indexes if others not in missing_pcgs]

    picked_fasta = path.join(basedir, f'{prefix}.picked.fa')
    SeqIO.write(picked_seq, picked_fasta, 'fasta')
    print(f'pcgs found : {found_pcgs}\nMissing pcgs : {missing_pcgs}')
    return picked_fasta, found_pcgs, missing_pcgs

@profiling
def nhmmer_search(fasta_file=None, thread_number=None, nhmmer_profile=None,
                  prefix=None, basedir=None):

    # Call nhmmer
    hmm_out = os.path.join(basedir, f'{prefix}.nhmmer.out')
    hmm_tbl = os.path.join(basedir, f'{prefix}.nhmmer.tblout')
    shell_call('nhmmer', o=hmm_out, tblout=hmm_tbl,
               cpu=thread_number, appending=[nhmmer_profile, fasta_file])

    # Process data to pandas readable table
    hmm_tbl_pd = f'{hmm_tbl}.readable'
    with open(hmm_tbl, 'r') as fin, open(hmm_tbl_pd, 'w') as fout:
        for line in fin:
            striped = line.strip()
            splitted = striped.split()
            # Dispose the description of genes, god damned nhmmer...
            print(' '.join(splitted[:15]), file=fout)

    # Read table with pandas
    hmm_frame = pandas.read_table(hmm_tbl_pd, comment='#', delimiter=' ',
                                  names=[
                                      'target', 'accession1', 'query',
                                      'accession2', 'hmmfrom', 'hmm to',
                                      'alifrom', 'alito', 'envfrom', 'envto',
                                      'sqlen', 'strand', 'e', 'score',
                                      'bias'
                                  ])
    hmm_frame = hmm_frame.drop(columns=['accession1', 'accession2'])

    # Deduplicate multiple hits on the same gene of same sequence
    hmm_frame = hmm_frame.drop_duplicates(
        subset=['target', 'query'], keep='first')
    hmm_frame.to_csv(f'{hmm_tbl}.dedup.csv', index=False)

    return hmm_frame


def filter_taxanomy(taxa=None, fasta_file=None, hmm_frame: pandas.DataFrame = None, basedir=None,
                    prefix=None, dbfile=None, gene_code=9, relaxing=0):

    # Extract sequences from input fasta file according to hmm frame

    seqs = [record
            for record in SeqIO.parse(fasta_file, 'fasta')
            if record.id in set(hmm_frame['target'])
            ]
    if not seqs:
        raise Exception("Parsed fasta file is empty!")

    hmm_fa = path.join(basedir, f'{prefix}.hmm.filtered.fa')
    with open(hmm_fa, 'w') as f:
        SeqIO.write(seqs, f, 'fasta')

    # Do tblastn to search out the possible taxanomy of the gene
    blast_file = tk.tblastn(dbfile=dbfile, infile=hmm_fa,
                            genetic_code=gene_code, basedir=basedir, prefix=prefix)
    blast_frame, _ = tk.blast_to_csv(blast_file)
    blast_frame = tk.wash_blast_results(blast_frame)

    # Drop the sequences which don't have even a gene related to taxa
    by_seqid = dict(tuple(blast_frame.groupby(['sseq'])))
    to_drop = []
    for key, frame in by_seqid.items():
        is_in = False
        for _, row in frame.iterrows():
            qseq = str(row.qseq).split('_')
            taxa_name = ' '.join([qseq[4], qseq[5]])
            taxa_rank = get_rank(taxa_name)
            required_rank = get_rank(taxa)
            required_id = ncbi.get_name_translator([taxa])[taxa][0]
            required_class = ncbi.get_rank([required_id])[required_id]
            required_index = rank_list.index(required_class)
            # Get last index for the matching rank
            matched_rank = max(idx
                               for idx, (tax_id, tax_name), (required_id, required_name)
                               in enumerate(zip(taxa_rank, required_rank))
                               if required_name == taxa_name != 'NA')
            if matched_rank + relaxing >= required_index:
                is_in = True
                break
        if not is_in:
            to_drop.append(key)

    filtered_frame = hmm_frame[~hmm_frame['target'].isin(to_drop)]
    return filtered_frame
