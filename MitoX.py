"""
MitoX.py
========

Copyright (c) 2019-2020 Li Junyu <2018301050@szu.edu.cn>.

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
    from utility.parser import *
    from arguments import *
except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoX installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

ncbi = NCBITaxa()

# Command processing
desc = """
Description

    MitoX - A rewrite toolkit of its ancestor MitoZ for faster and better 
    mitochondrial assembly, annotation and visualization, and for expandability and more.

Version
    1.0
"""

@parse_func(func_help='filter out unqualified reads from fastq',
            parents=[fastq_parser, filter_parser, universal_regulator])
@arg_prop(dest='seq_size', help='how many sequences will be filtered out.', arg_type=int)
def filter(args):

    dest = args.temp_folder if args.__calling == 'filter' else args.result_folder
    cleanq1 = path.join(dest, 'filtered.1.fq.gz')
    if args.fastq2 is not None:
        cleanq2 = path.join(dest, 'filtered.2.fq.gz')

    filtered1 = filtered2 = None

    from filter.filter import filter_pe, filter_se

    if args.fastq2 is None:
        filtered1 = filter_se(fqiabs=f'"{args.fastq1}"', fqoabs=f'"{args.cleanq1}"', Ns=args.Ns_valve,
                              quality=args.quality_valve, limit=args.percentage_valve, start=args.start,
                              end=args.end, seq_size=args.seq_size)
    else:
        filtered1, filtered2 = filter_pe(fq1=f'"{args.fastq1}"', fq2=f'"{args.fastq2}"',
                                         o1=f'"{cleanq1}"', o2=f'"{cleanq2}"',
                                         a1=f'"{args.adapter1}"' if args.adapter1 is not None else None,
                                         a2=f'"{args.adapter2}"' if args.adapter2 is not None else None,
                                         dedup=args.deduplication,
                                         mis=args.adapter_mismatch, ali=args.adapter_length, start=args.start,
                                         end=args.end, n=args.Ns_valve, q=args.quality_valve, l=args.percentage_valve,
                                         seq_size=args.seq_size)
    return filtered1, filtered2

@parse_func(func_help='assemble from input fastq reads, output contigs',
            parents=[universal_parser, fastq_parser, assembly_parser])
def assemble(args):

    from assemble.assemble import assemble

    assembled_contigs = assemble(fastq1=args.fastq1, fastq2=args.fastq2, result_dir=args.result_folder,
                                 temp_dir=args.temp_folder, work_prefix=args.workname, uselist=args.use_list,
                                 kmin=args.kmer_min, kmax=args.kmer_max, kstep=args.kmer_step, klist=args.kmer_list,
                                 no_mercy=args.no_mercy, disable_acc=args.disable_acc, prune_level=args.prune_level,
                                 prune_depth=args.prune_depth, clean_temp=args.clean_temp, threads=args.threads)

    return assembled_contigs


@parse_func(func_help='search for the most possible mitochondrial sequences from assembled data',
            parents=[universal_parser, fasta_parser, fastq_parser, search_parser, saa_parser])
def findmitoscaf(args):
    # TODO:To fill the blanks of findmitoscaf method
    pass


@parse_func(func_help='annotate PCGs, tRNA and rRNA genes',
            parents=[universal_parser, fasta_parser, annotation_parser, saa_parser, fastq_parser])
@arg_prop(dest='depth_file', help=argparse.SUPPRESS)
@arg_prop(dest='topology', choices=['linear', 'circular'], help='if the sequences are circular')
def annotate(args):
    # TODO:To fill the blanks of annotate method
    pass


@parse_func(func_help='visualization of GenBank file')
def visualize(args):
    # TODO:To fill the blanks of visualize method
    pass


@parse_func(func_help='run all the methods',
            parents=[universal_parser, assembly_parser, fastq_parser, filter_parser,
                     search_parser, saa_parser, annotation_parser])
@arg_prop(dest='disable_filter', help='filter will be not enabled if this switched on')
@arg_prop(dest='topology', choices=['linear', 'circular'], help=argparse.SUPPRESS)
@arg_prop(dest='seq_size', help='how many sequences will be filtered out.', arg_type=int)
def all(args):

    # Go filtering
    if not args.disable_filter:
        args.fastq1, args.fastq2 = filter(args=args)

    contigs_file = assemble(args)

    # TODO:Finish findmitoscaf methods.


# Entry starts at here
if __name__ == '__main__':
    parser = freeze_arguments('MitoX', desc)
    parse_then_call(parser)
