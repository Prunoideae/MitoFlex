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

import pandas
from Bio import SeqIO
from ete3 import NCBITaxa

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call
    from utility.profiler import profiling
    from utility.seq import compile_seq, decompile
    from annotation.annotation_tookit import *
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

ncbi = NCBITaxa()


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

    return rank_dict


@profiling
def findmitoscaf():
    pass


@profiling
def nhmmer_search(fasta_file=None, thread_number=None, nhmmer_profile=None,
                  prefix=None, basedir=None, gene_code=9, dbfile=None, taxa=None):

    # Decompress zipped files
    input_file = fasta_file
    if input_file.endswith(".gz"):
        input_file = f'{prefix}.tmp.nhmmer.dat'
        shell_call('gzip', '-dc', fasta_file, '>', input_file)

    # Call nhmmer
    hmm_out = os.path.join(basedir, f'{prefix}.nhmmer.out')
    hmm_tbl = os.path.join(basedir, f'{prefix}.nhmmer.tblout')
    shell_call('nhmmer', o=hmm_out, tblout=hmm_tbl,
               cpu=thread_number, appending=[nhmmer_profile, input_file])

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
    hmmframe = hmm_frame.drop_duplicates(
        subset=['target', 'query'], keep='first')
    hmm_frame.to_csv(f'{hmm_tbl}.dedup.csv')

    return hmm_frame


def filter_taxanomy(taxa=None, fasta_file=None, hmm_frame=None, basedir=None,
                    prefix=None, dbfile=None, gene_code=9):

    # Extract sequences from input fasta file according to hmm frame
    records = SeqIO.parse(fasta_file, 'fasta')
    seqs = [record for record in records if record.name in set(
        hmm_frame['target'])]
    if not seqs:
        return None, None

    hmm_fa = path.join(basedir, f'{prefix}.hmm.filtered.fa')
    with open(hmm_fa, 'w') as f:
        SeqIO.write(seqs, f, 'fasta')

    blast_file = tblastn(dbfile=dbfile, infile=hmm_fa,
                         genetic_code=gene_code, basedir=basedir, prefix=prefix)
    blast_frame, _ = blast_to_csv(blast_file)
    blast_frame = wash_blast_results(blast_frame)
    wise_frame, queries, database = genewise(basedir=basedir, prefix=prefix,
                                             wises=blast_frame, infile=hmm_fa, dbfile=dbfile)
    final_file = path.join(basedir, f'{prefix}.genewise.cds.fa')
    if not collect_result(final_file, wise_frame, queries, database):
        return None, None

    required_rank = get_rank(taxa)
    if not required_rank:
        return None, None

    valids = []
    records = SeqIO.parse(final_file, 'fasta')
    for record in records:
        if record.description.startswith(record.id):
            record.description.replace(record.id + " ", '', 1)
        traits = decompile(record.description)
        qseq = traits['qseq'].split('_')
        sseq = traits['sseq']
        qspecies = ' '.join(qseq[4:6])
        if qspecies.endswith('.'):
            qspecies = qspecies.split('.')[0]

        ranks = get_rank(qspecies)
        valids.append((qspecies, sseq, ranks))

    valids.append((None, None, required_rank))
    valids = list(set(valids))
    valids.sort(key=lambda x: '|'.join(str(i[1] for i in x[3].items())))
