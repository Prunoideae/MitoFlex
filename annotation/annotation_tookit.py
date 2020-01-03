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

    blast_frame = pandas.read_table(
        out_blast, delimiter='\t',
        names=['qseq', 'sseq', 'ident',
               'length', 'mismatch', 'gap',
               'qstart', 'qend', 'sstart',
               'send', 'evalue', 'score'])

    # Delete duplicated and unqualified results, then add extra informations about results.
    blast_frame = blast_frame[blast_frame.ident < ident]
    blast_frame = blast_frame.drop_duplicates(keep='first')
    blast_frame = blast_frame[blast_frame.score > 25]
    blast_frame['qmax'] = blast_frame.groupby('qseq')['qend'].transform(
        lambda x: max(x) if x.count() > 2 else x)
    blast_frame = blast_frame[blast_frame.qstart - blast_frame.qend
                              >= blast_frame.qmax*0.25]
    blast_frame['qseq'] = blast_frame.groupby('qseq')['qseq'].transform(
        lambda s: [f'{x}_{idx}' for idx, x in enumerate(s)]
    )

    # For logging purpose
    out_blast_csv = f'{prefix}.blast.csv'
    blast_frame.to_csv(out_blast_csv)

    return out_blast_csv, blast_frame


def genewise(basedir=None, prefix=None, blast_frame: pandas.DataFrame = None):

    # Filter out the most important sequences for genewiseF
    blast_frame['plus'] = blast_frame.send - blast_frame.sstart > 0
    blast_frame['sstart'], blast_frame['send'] = np.where(
        blast_frame['sstart'] > blast_frame['send'],
        [blast_frame['send'], blast_frame['sstart']],
        [blast_frame['sstart'], blast_frame['send']]
    )

    by_sseq = dict(tuple(blast_frame.groupby(['sseq', 'plus'])))
    by_sseq = {key: value.sort_values('sstart')
               for key, value in by_sseq.items()}

    cluster = {}

    for identifier, frame in by_sseq.items:
        cluster = {}
        for row in frame.iterrows():
            pass
