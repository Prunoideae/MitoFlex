"""
annotation_toolkit.py
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
from os import path
import sys

import pandas
import numpy as np
from Bio import SeqIO

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import concat_command, direct_call, shell_call
    from utility.profiler import profiling
except Exception as iden:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)

# Truncates all the -- to - to suit blast's parsing style


def truncated_call(*args, **kwargs):
    direct_call(concat_command(*args, **kwargs).replace('--', '-'))


def tblastn(dbfile=None, infile=None, genetic_code=9, basedir=None,
            prefix=None, ident=30):
    # Since 'in' is the preserve keyword
    truncated_call('makeblastdb', '-in', dbfile, dbtype='nucl')

    # Go tblastn
    out_blast = path.join(path.abspath(basedir), f'{prefix}.blast')
    truncated_call('tblastn', useconv=False, evalue='1e-5', outfmt=6,
                   seg='no', db_gencode=genetic_code, db=infile,
                   query=dbfile)

    return out_blast


def blast_to_csv(blast_file, ident=30, score=25):
    blast_frame = pandas.read_table(
        blast_file, delimiter='\t',
        names=['qseq', 'sseq', 'ident',
               'length', 'mismatch', 'gap',
               'qstart', 'qend', 'sstart',
               'send', 'evalue', 'score'])

    # Delete duplicated and unqualified results, then add extra informations about results.
    blast_frame = blast_frame.drop_duplicates(keep='first')
    blast_frame = blast_frame[blast_frame.ident > ident]
    blast_frame = blast_frame[blast_frame.score > score]
    blast_frame['qmax'] = blast_frame.groupby('qseq')['qend'].transform(
        lambda x: max(x) if x.count() > 2 else x)
    blast_frame = blast_frame[blast_frame.qend - blast_frame.qstart
                              >= blast_frame.qmax*0.25]

    # For logging purpose
    out_blast_csv = f'{blast_file}.csv'
    blast_frame.to_csv(out_blast_csv)

    return blast_frame, out_blast_csv


# Filter out the most important sequences for genewise
def wash_blast_results(blast_frame: pandas.DataFrame = None, ident=0.5):
    blast_frame['plus'] = blast_frame.send - blast_frame.sstart > 0
    blast_frame['sstart'], blast_frame['send'] = np.where(
        blast_frame['sstart'] > blast_frame['send'],
        [blast_frame['send'], blast_frame['sstart']],
        [blast_frame['sstart'], blast_frame['send']]
    )

    by_sseq = dict(tuple(blast_frame.groupby(['sseq', 'plus'])))
    by_sseq = {key: value.sort_values('sstart')
               for key, value in by_sseq.items()}

    results = []

    for identifier, frame in by_sseq.items():
        # Find all highest score results which does not overlap with
        # any other sequences.
        while not frame.empty:
            highest = frame[frame.score == frame.score.max()].head(1)
            results.append(highest)

            max_len = int(highest.send - highest.sstart)+1
            max_start = int(highest.sstart)
            max_end = int(highest.send)

            frame = frame.drop(highest.index)
            cutoffs = np.min(max_len, frame.send - frame.sstart) * ident
            overlays = np.min(frame.send - max_start, max_end - frame.send)
            frame = frame.drop(overlays > cutoffs)

    wises = pandas.DataFrame(pandas.concat(results))
    return wises


def genewise(basedir=None, prefix=None, codon_table=None,
             blast_frame: pandas.DataFrame = None, infile=None,
             dbfile=None, ident=0.5):

    wises = wash_blast_results(blast_frame, ident)
    wises.to_csv(path.join(basedir, f'{prefix}.targets'))

    wisedir = path.join(basedir, 'genewise')
    dbdir = path.join(wisedir, 'sequences')
    query_dir = path.join(wisedir, 'queries')

    queries = {record.id: record
               for record in SeqIO.parse(infile, 'fasta')
               if record.id in set(wises.sseq)}

    for record in SeqIO.parse(dbfile, 'fasta'):
        if record.id in set(wises.qseq):
            SeqIO.write(record, path.join(dbdir, f'{record.id}.fa'), 'fasta')

    for index, wise in wises.iterrows():
        query_prefix = f'{wise.qseq}_{wise.sseq}_{wise.sstart}_{wise.send}'
        query_file = f'{query_prefix}.fa'
        query_result = f'{query_prefix}.genewise'

        SeqIO.write(queries[wise.sseq]
                    [wise.sstart-1:wise.send], query_file, 'fasta')

        truncated_call('genewise', condon=codon_table,
                       trev=not wise.plus, genesf=True,
                       gff=True, gum=True,
                       appending=[
                           path.join(dbdir, f'{wise.qseq}.fa'),
                           query_file,
                           '>', query_result
                       ])
