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
import subprocess

import pandas
import numpy as np
from Bio import SeqIO

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import concat_command, direct_call, shell_call
    from utility.profiler import profiling
    from utility.seq import compile_seq, decompile
except Exception as iden:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)

# Truncates all the -- to - to suit blast's parsing style


def truncated_call(*args, **kwargs):
    return direct_call(concat_command(*args, **kwargs).replace('--', '-'))


def tblastn(dbfile=None, infile=None, genetic_code=9, basedir=None,
            prefix=None, ident=30):

    # Make sure input path IS absolute path
    infile = path.abspath(infile)
    dbfile = path.abspath(dbfile)

    # Since 'in' is the preserve keyword
    truncated_call('makeblastdb', '-in', infile, dbtype='nucl')

    # Go tblastn
    out_blast = path.join(path.abspath(basedir), f'{prefix}.blast')
    truncated_call('tblastn', useconv=False, evalue='1e-5', outfmt=6,
                   seg='no', db_gencode=genetic_code, db=infile,
                   query=dbfile, appending=['>', out_blast])

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
    blast_frame.to_csv(out_blast_csv, index=False)

    return blast_frame, out_blast_csv


# Filter out the most important sequences for genewise
def wash_blast_results(blast_frame: pandas.DataFrame = None, cutoff=0.5):
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

    for frame in by_sseq.values():
        # Find all highest score results which does not overlap with
        # any other sequences.
        while not frame.empty:
            highest = frame[frame.score == frame.score.max()].head(1)
            results.append(highest)

            max_len = int(highest.send - highest.sstart)+1
            max_start = int(highest.sstart)+1
            max_end = int(highest.send)

            frame = frame.drop(highest.index)
            cutoffs = np.minimum(max_len, frame.send - frame.sstart) * cutoff
            overlays = np.minimum(frame.send - max_start,
                                  max_end - frame.sstart)
            frame = frame[overlays <= cutoffs]

    return pandas.concat(results)


def genewise(basedir=None, prefix=None, codon_table=None,
             wises: pandas.DataFrame = None, infile=None,
             dbfile=None, cutoff=0.5):

    if codon_table is None:
        codon_table = path.join(bin_dir, 'codon_InverMito.table')

    wisedir = path.join(basedir, 'genewise')
    dbdir = path.join(wisedir, 'sequences')
    query_dir = path.join(wisedir, 'queries')

    queries = {record.id: record
               for record in SeqIO.parse(infile, 'fasta')
               if record.id in set(wises.sseq)}

    dbparsed = {record.id: record
                for record in SeqIO.parse(dbfile, 'fasta')
                if record.id in set(wises.qseq)}

    for idx, record in dbparsed.items():
        SeqIO.write(record, path.join(dbdir, f'{idx}.fa'), 'fasta')

    wises = wises.assign(
        wise_min_start=np.nan, wise_max_end=np.nan,
        wise_shift=np.nan, wise_cover=np.nan)

    wise_cfg_dir = path.join(bin_dir, 'wisecfg')
    env_var = dict(os.environ)
    env_var["WISECONFIGDIR"] = wise_cfg_dir

    for _, wise in wises.iterrows():
        query_prefix = f'{wise.qseq}_{wise.sseq}_{wise.sstart}_{wise.send}'
        query_file = path.join(query_dir, f'{query_prefix}.fa')

        SeqIO.write(queries[wise.sseq]
                    [wise.sstart-1:wise.send], query_file, 'fasta')

        result = subprocess.check_output(
            concat_command('genewise', codon=codon_table,
                           trev=not wise.plus, genesf=True, gff=True, gum=True,
                           appending=[
                               path.join(dbdir, f'{wise.qseq}.fa'), query_file]
                           ).replace("--", '-'), env=env_var, shell=True).decode('utf-8')

        # Parse the results
        splited = result.split('//\n')
        info = splited[0].split('\n')[1].split()
        wise.wises.at[idx, 'wise_cover'] = (
            int(info[3]) - int(info[2]) + 1)/len(dbparsed[wise.qseq])
        wise_result = [x.split('\t')
                       for x in splited[2].split('\n')[:-1]
                       if x.split('\t')[2] == 'cds']
        for x in wise_result:
            # Fix the actual position of seq
            x[3] = str(int(x[3]) + wise.sstart-1)
            x[4] = str(int(x[4]) + wise.sstart-1)
        wise_result.sort(key=lambda x: int(x[3]))
        wises.at[idx, 'wise_shift'] = sum(
            x[2] == 'match' for x in wise_result)-1
        wises.at[idx, 'wise_min_start'] = min(int(x[3]) for x in wise_result)
        wises.at[idx, 'wise_max_end'] = max(int(x[4]) for x in wise_result)

    wises.to_csv(path.join(basedir, f'{prefix}.wise.csv'), index=False)
    return wises, queries, dbparsed


def collect_result(output_file=None, wises: pandas.DataFrame = None, queries=None, data=None):
    result_seq = []
    for _, wise in wises.iterrows():
        trait_string = compile_seq({
            'qseq': wise.qseq,
            'sseq': wise.sseq,
            'qstart': wise.qstart,
            'qend': wise.qend,
            'sstart': wise.sstart,
            'send': wise.send,
            'plus': wise.plus
        })

        seq = queries[wise.sseq][wise.sstart-1:wise.send]
        seq.description = trait_string
        if wise.plus is 0:
            seq.seq = seq.seq.reverse_complement()
        result_seq.append(seq)
    return SeqIO.write(result_seq, output_file, 'fasta')


def reloc_genes(fasta_file=None, wises: pandas.DataFrame = None, code=9):
    wise_seqs = {x.id: x for x in SeqIO.parse(fasta_file, 'fasta')}
    wises.assign(start_real=np.nan, end_real=np.nan)
    for _, wise in wises.iterrows():
        start_real = end_real = -1
        seq = wise_seqs[wise.sseq][wise.sstart-29:wise.send+30]
        if not wise.plus:
            seq = seq.reverse_complement()
        try:
            trans = seq.translate(9, cds=True)
        except:
            trans = seq.translate(9)
        # Finding stop
        if trans.find('*') != -1:
            offset = trans.find('*') * 3
            end_real = (wise.sstart+1 + offset
                        if wise.plus else
                        wise.send - offset)
        else:
            mercy = seq[-30:] if wise.plus else seq[:30].reverse_complement()
            offset = mercy.find('TA')
            if offset == -1:
                offset = mercy.find('T')

            if offset != -1:
                end_real = (
                    wise.send+30 + offset if wise.plus else wise.sstart - 28 - offset)

        # Finding start
        # Find the last M of mercied part of translated seq, where it should be the start codon
        mercy = trans[:10:-1].find('M')
        if mercy != -1:
            offset = (10-mercy) * 3
            start_real = wise.sstart + offset - 29 if wise.plus else wise.send + 30 + offset

        wise.start_real, wise.end_real = (
            start_real, end_real) if wise.plus else (end_real, start_real)

    return wises
