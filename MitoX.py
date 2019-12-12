"""
MitoX.py
========

Copyright (c) 2019-2020 Henry Lee <2018301050@szu.edu.cn>.

This file is part of MitoX.

MitoX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoX is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoX.  If not, see <http://www.gnu.org/licenses/>.

"""

import argparse
import sys
import os
from os import path
import re
import subprocess
import time
from glob import glob

try:
    import Bio
    from Bio import SeqIO
    from Bio import SeqRecord
    from ete3 import NCBITaxa
    from helper import *
except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoX installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

ncbi = NCBITaxa()

# Environment intialization
python3 = sys.executable
work_dir = os.getcwd()

filter_dir = os.path.join(sys.path[0], 'filter')
assemble_dir = os.path.join(sys.path[0], 'assemble')
profile_dir = os.path.join(sys.path[0], 'profiles')
visualize_dir = os.path.join(sys.path[0], 'visualize')
findmitoscaf_dir = os.path.join(sys.path[0], 'findmitoscaf')
misc_dir = os.path.join(sys.path[0], 'misc')
annotate_dir = os.path.join(sys.path[0], 'annotate')

# Parameters processing

# Universal arguments
universal_parser, universal_group = register_group('Universal arguments', [
    {
        'name': 'workname',
        'required': True,
        'help': 'work output name.'
    },
    {
        'name': 'threads',
        'default': '8',
        'help': 'thread numbers.'
    },
    {
        'name': 'clean-temp',
        'default': False,
        'help': 'remove temporal files and folders after work done.'
    },
    {
        'name': 'basedir',
        'default': None,
        'help': 'result and temp folder will be generated in this instead of current directory.'
    }
])

# Fastq arguments
fastq_parser, fastq_group = register_group('Fastq arguments', [
    {
        'name': 'fastq1',
        'default': None,
        'meta': 'file',
        'help': 'fastq file 1 input.'
    },
    {
        'name': 'fastq2',
        'default': None,
        'meta': 'file',
        'help': 'fastq file 2 input.'
    },
    {
        'name': 'fastq-alter-format',
        'default': False,
        'action': 'store_true',
        'help': 'using FASTQ-like (Q+64) format or Q+33 file format.'
    },
    {
        'name': 'fastq-read-length',
        'default': 150,
        'help': '''read length of fastq reads, used in MEGAHIT and bwa,
        it must be at least 71 bp.'''
    }
])

# Fasta group
fasta_parser, fasta_group = register_group('Fasta arugments', [
    {
        'name': 'fastafile',
        'meta': 'file',
        'help': 'fasta file'
    }
])

# Filter arguments
filter_parser, filter_group = register_group('Filter argumetns', [
    {
        'name': 'adapter1',
        'meta': 'file',
        'help': 'input adapter list file 1 to filter adapter contamination.'
    },
    {
        'name': 'adapter2',
        'meta': 'file',
        'help': 'input adapter list file 1 to filter adapter contamination.'
    },
    {
        'name': 'deduplication',
        'default': False,
        'action': 'store_true',
        'help': 'fitler duplication caused by adapter ligation if switched on.'
    },
    {
        'name': 'adapter-mismatch',
        'default': 3,
        'help': 'cut-off adapter mismatch bases.'
    },
    {
        'name': 'adapter-length',
        'default': 15,
        'help': 'cut-off adapter align length.'
    },
    {
        'name': 'Ns-valve',
        'default': 10,
        'help': 'the read will be discarded if more than X Ns in the sequence.'
    },
    {
        'name': 'quality-valve',
        'default': 55,
        'help': 'bases with quality under this will be counted as a bad base.'
    },
    {
        'name': 'percentage-valve',
        'default': 0.2,
        'help': 'a read have a percentage of bad bases over the total length bigger this will be discarded.'
    },
    {
        'name': 'keep-region',
        'default': '-1,-1',
        'meta': 'beg,end',
        'help': 'only the BEG and the END will be read, leave blank for full length.'
    }
])

# Assembly group
assembly_parser, assembly_group = register_group('Assembly arguments', [
    {
        'name': 'insert-size',
        'default': 150,
        'help': 'insert size of input fastq files.'
    },
    {
        'name': 'use-list',
        'default': False,
        'action': 'store_true',
        'help': 'a k-mer list will be used if this switched on.'
    },
    {
        'name': 'kmer-list',
        'default': '21,29,39,59,79,99,119,141',
        'help': 'list of kmer to use in sDBG building, all length must be odd.'
    },
    {
        'name': 'kmer-min',
        'default': 21,
        'help': 'the minimum length kmer used in MEGAHIT\'s multiple kmer assemble strategy.'
    },
    {
        'name': 'kmer-max',
        'default': 141,
        'help': 'the maximum length kmer used in MEGAHIT\'s multiple kmer assemble strategy.'
    },
    {
        'name': 'kmer-step',
        'default': 12,
        'help': 'increment of kmer size of each iteration (<= 28), must be even number.'
    },
    {
        'name': 'no-mercy',
        'default': False,
        'help': 'mercy edges are NOT allowed in sDBG if switched on.',
        'action': 'store_true'
    },
    {
        'name': 'prune-level',
        'default': 2,
        'choices': list(range(0, 4)),
        'help': 'strength of low depth pruning.'
    },
    {
        'name': 'prune-depth',
        'default': 2,
        'help': 'remove unitigs with avg kmer depth less than this value.'
    },
    {
        'name': 'disable-acc',
        'default': False,
        'action': 'store_true',
        'help': 'force disable the hardware accerlation if there\'s some problem in using it, only use as a last resort.'
    }
])

# Search mitochondrial gene group
search_parser, search_group = register_group('Search mitochondrial sequences arguments', [
    {
        'name': 'filter-taxa',
        'default': False,
        'action': 'store_true',
        'help': 'turn on to filter out sequences not match the clade given.'
    },
    {
        'name': 'min-abundance',
        'default': 10,
        'type': float,
        'help': 'the minimum abundance of the sequence required in filtering.',
    },
    {
        'name': 'required-taxa',
        'default': 'Platyhelminthes',
        'help': 'taxa sequences other than this will be filtered out.'
    },
    {
        'name': 'taxa-tolerance',
        'default': 0,
        'choices': list(range(7)),  # 0-6
        'help': 'a tolerance for non-target sequences.'
    }
])

# Search and annotation mitochondrial gene group
saa_parser, saa_group = register_group('Search and annotate arguments', [
    {
        'name': 'genetic-code',
        'default': 9,
        'help': 'genetic code table to be used in the run.'
    },
    {
        'name': 'clade',
        'default': 'Platyhelminthes-flatworms',
        'choices': ['Chordata', "Platyhelminthes-flatworms"],
        'help':'which clade\'s nhmmer profile and cds will be used in the run.'
    }
])

# Annotation group
annotation_parser, annotation_group = register_group('Annotation arguments', [
    {
        'name': 'disable-annotation',
        'action': 'store_false',
        'default': True,
        'help': 'switching this on will disable annotation.'
    },
    {
        'name': 'species-name',
        'default': 'Test sp.',
        'help': 'species name to use in genbank file.'
    }
])


# Command processing
desc = """
Description

    MitoX - A rewrite toolkit of its ancestor MitoZ for faster and better 
    mitochondrial assembly, annotation and visualization, and for expandability and more.

Version
    1.0
"""

@parse_func(func_help='filter out unqualified reads from fastq',
            parents=[fastq_parser, filter_parser])
@arg_prop(dest='workname', help='in where the filtered data will be put in.', required=True)
@arg_prop(dest='basedir', help='result will be in this directory else than the current')
@arg_prop(dest='seq_size', help='how many sequences will be filtered out.', arg_type=int)
def filter(workname=None, seq_size=None, basedir=None,

           fastq1=None, fastq2=None, fastq_alter_format=False, fastq_read_length=150,

           adapter1=None, adapter2=None, deduplication=False,
           adapter_mismatch=3, adapter_length=15, keep_region='-1,-1', Ns_valve=10, quality_valve=50,
           percentage_valve=0.2
           ):

    if fastq1 is None and fastq2 is not None:
        fastq1, fastq2 = fastq2, fastq1

    elif fastq1 is None and fastq2 is None:
        collected_args['filter']['parser'].print_help()
        sys.exit('Both fastq input file is not specified. Exiting.')

    work_folder = path.abspath(
        path.join(work_dir if basedir is None else basedir, workname))
    temp_folder = path.join(work_folder, workname + '.tmp')

    cleanq1 = path.join(temp_folder, 'filtered.1.fq.gz')
    cleanq2 = path.join(temp_folder, 'filtered.2.fq.gz')

    # Validates the directories
    try:
        os.makedirs(path.dirname(cleanq1), exist_ok=True)
        if fastq2 is not None:
            os.makedirs(path.dirname(cleanq2), exist_ok=True)
    except Exception as identifier:
        sys.exit(
            'Error occured when validating the directories, please check your permissions or things could be realted.')

    fastq1 = path.abspath(fastq1)
    if fastq2 is not None:
        fastq2 = path.abspath(fastq2)

    filtered1 = None
    filtered2 = None

    start, end, *_ = [int(x) if int(x) > -1 else None
                      for x in keep_region.split(',')]

    from filter.filter import filter_pe, filter_se

    if fastq1 is not None and fastq2 is None:
        filtered1 = filter_se(fqiabs=f'"{fastq1}"', fqoabs=f'"{cleanq1}"', Ns=Ns_valve,
                              quality=quality_valve, limit=percentage_valve, start=start,
                              end=end, seq_size=seq_size)
    else:
        filtered1, filtered2 = filter_pe(fq1=f'"{fastq1}"', fq2=f'"{fastq2}"',
                                         o1=f'"{cleanq1}"', o2=f'"{cleanq2}"',
                                         a1=f'"{adapter1}"' if adapter1 is not None else None,
                                         a2=f'"{adapter2}"' if adapter2 is not None else None,
                                         dedup=deduplication,
                                         mis=adapter_mismatch, ali=adapter_length, start=start,
                                         end=end, n=Ns_valve, q=quality_valve, l=percentage_valve,
                                         seq_size=seq_size)
    return filtered1, filtered2

@parse_func(func_help='assemble from input fastq reads, output contigs',
            parents=[universal_parser, fastq_parser, assembly_parser])
def assemble(workname=None, threads=8, clean_temp=False, basedir=None,

             fastq1=None, fastq2=None, fastq_alter_format=None, fastq_read_length=125,

             insert_size=150, kmer_min=21, kmer_max=141, kmer_step=12, kmer_list=None, use_list=False, no_mercy=False,
             prune_level=2, prune_depth=2, disable_acc=False,
             ):

    if fastq1 is None and fastq2 is not None:
        fastq1, fastq2 = fastq2, fastq1

    if not fastq1 and not fastq2:
        sys.exit('Both fastq1 and fastq2 are not specified!')

    if (use_list):
        if 0 in [int(x) % 2 for x in kmer_list.split(',')]:
            sys.exit('K-mer specified should always be odd numbers!')
    else:
        if 0 in [int(x) % 2 for x in range(kmer_min, kmer_max, kmer_step)]:
            sys.exit('K-mer specified should always be odd numbers!')

    if (threads <= 0):
        sys.exit('Unreasonable thread count specified!')

    work_folder = path.abspath(
        path.join(work_dir if basedir is None else basedir, workname))
    result_folder = path.join(work_folder, workname + '.result')
    temp_folder = path.join(work_folder, workname + '.tmp')

    try:
        os.makedirs(result_folder, exist_ok=True)
        os.makedirs(temp_folder, exist_ok=True)
    except Exception as identifier:
        print(identifier)
        sys.exit(
            'Error occured when validating the directories, please check your permissions or things could be realted.')

    from assemble.assemble import assemble

    assembled_contigs = assemble(fastq1=fastq1, fastq2=fastq2, result_dir=result_folder,
                                 temp_dir=temp_folder, work_prefix=workname, uselist=use_list,
                                 kmin=kmer_min, kmax=kmer_max, kstep=kmer_step, klist=kmer_list,
                                 no_mercy=no_mercy, disable_acc=disable_acc, prune_level=prune_level,
                                 prune_depth=prune_depth, clean_temp=clean_temp, threads=threads)

    return assembled_contigs


@parse_func(func_help='search for the most possible mitochondrial sequences from assembled data',
            parents=[universal_parser, fasta_parser, fastq_parser, search_parser, saa_parser])
def findmitoscaf(workname=None, threads=8, clean_temp=False,

                 fastq1=None, fastq2=None, fastq_alter_format=None, fastq_read_length=150,

                 filter_taxa=False, min_abundance=10, required_taxa='Platyhelminthes', taxa_tolerance=0,

                 genetic_code=9, clade='Platyhelminthes-flatworms',

                 fastafile=None):
    # TODO:To fill the blanks of findmitoscaf method
    pass


@parse_func(func_help='annotate PCGs, tRNA and rRNA genes',
            parents=[universal_parser, fasta_parser, annotation_parser, saa_parser])
@arg_prop(dest='fastq1', help='fastq file 1 to draw depth')
@arg_prop(dest='fastq2', help='fastq file 2 to draw depth')
@arg_prop(dest='depth_file', help=argparse.SUPPRESS)
@arg_prop(dest='topology', choices=['linear', 'circular'], help='if the sequences are circular')
def annotate(fastq1=None, fastq2=None, depth_file=None, topology='linear',

             workname=None, threads=8, clean_temp=False,

             fastafile=None,

             disable_annotation=False, species_name='Test sp.',

             genetic_code=9, clade='Platyhelminthes-flatworms'):
    # TODO:To fill the blanks of annotate method
    pass


@parse_func(func_help='visualization of GenBank file')
def visualize():
    # TODO:To fill the blanks of visualize method
    pass


@parse_func(func_help='run all the methods',
            parents=[universal_parser, assembly_parser, fastq_parser, filter_parser,
                     search_parser, saa_parser, annotation_parser])
@arg_prop(dest='disable_filter', help='filter will be not enabled if this switched on')
@arg_prop(dest='topology', choices=['linear', 'circular'], help=argparse.SUPPRESS)
@arg_prop(dest='seq_size', help='how many sequences will be filtered out.', arg_type=int)
def all(topology='linear', disable_filter=False, seq_size=None,

        workname=None, threads=8, clean_temp=False, basedir=None,

        fastq1=None, fastq2=None, fastq_alter_format=None, fastq_read_length=150,

        adapter1=None, adapter2=None, deduplication=False,
        adapter_mismatch=3, adapter_length=15, keep_region=None, Ns_valve=10, quality_valve=50,
        percentage_valve=20,

        insert_size=150, kmer_min=31, kmer_max=63, kmer_step=12, kmer_list=None, use_list=False, no_mercy=False,
        prune_level=2, prune_depth=2, disable_acc=False,

        filter_taxa=False, min_abundance=10, required_taxa='Platyhelminthes', taxa_tolerance=0,

        genetic_code=9, clade='Platyhelminthes-flatworms',

        disable_annotation=False, species_name='Test sp.'
        ):

    dataq1, dataq2 = fastq1, fastq2
    working_folder = path.abspath(
        path.join(work_dir if basedir is None else basedir, workname))
    temp_folder = path.join(working_folder, workname + '.tmp')

    # Go filtering
    if not disable_filter:
        dataq1, dataq2 = filter(workname=workname, seq_size=seq_size, fastq1=fastq1, fastq2=fastq2, adapter1=adapter1,
                                adapter2=adapter2, deduplication=deduplication,
                                adapter_mismatch=adapter_mismatch, adapter_length=adapter_length,
                                keep_region=keep_region, Ns_valve=Ns_valve, quality_valve=quality_valve,
                                percentage_valve=percentage_valve, basedir=basedir)
    elif dataq1 is None and dataq2 is not None:
        dataq1, dataq2 = dataq2, dataq1

    contigs_file = assemble(workname=workname, threads=threads, clean_temp=clean_temp, basedir=basedir,
                            fastq1=dataq1, fastq2=dataq2, fastq_alter_format=fastq_alter_format,
                            fastq_read_length=fastq_read_length, insert_size=insert_size,
                            kmer_min=kmer_min, kmer_max=kmer_max, kmer_step=kmer_step, kmer_list=kmer_list,
                            use_list=use_list, no_mercy=no_mercy, prune_level=prune_level, prune_depth=prune_depth,
                            disable_acc=disable_acc)

    # TODO:Finish findmitoscaf methods.


# Entry starts at here
if __name__ == '__main__':
    parser = freeze_arguments('MitoX', desc)
    parse_then_call(parser)
