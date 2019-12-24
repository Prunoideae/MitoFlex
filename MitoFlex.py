#!/usr/bin/env python3

"""
MitoFlex.py
========

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

from glob import glob
import time
import subprocess
import re
from os import path
import os
import sys
import argparse

if sys.version_info[0] < 3:
    sys.exit('Python 3 must be installed in current environment! Please check if any of your environment setup(like conda environment) is deactivated or wrong!')

try:
    import Bio
    from Bio import SeqIO
    from Bio import SeqRecord
    from ete3 import NCBITaxa
    from utility.parser import *
    from utility.profiler import profiling
    from arguments import *
except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoFlex installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

# Environment initialization
ncbi = NCBITaxa()

# Command processing
desc = """
Description

    MitoFlex - A rewritten toolkit of its ancestor MitoZ for faster and better 
    mitochondrial assembly, annotation and visualization, for expandability and more.

Version
    0.0

Citation

"""

@parse_func(func_help='filter out unqualified reads from fastq',
            parents=[universal_parser, fastq_parser, filter_parser])
@arg_prop(dest='seq_size', help='how many sequences will be filtered out.', arg_type=int)
def filter(args):

    dest = args.result_folder if args.__calling == 'filter' else args.temp_folder
    if not path.isabs(args.cleanq1):
        args.cleanq1 = path.abspath(path.join(dest, args.cleanq1))

    if args.cleanq2 is not None and not path.isabs(args.cleanq2):
        args.cleanq2 = path.abspath(path.join(dest, args.cleanq2))

    filtered1 = filtered2 = None

    from filter.filter import filter_pe, filter_se

    if args.fastq2 is None:
        filtered1 = filter_se(fqiabs=f'"{args.fastq1}"', fqoabs=f'"{args.cleanq1}"', Ns=args.Ns_valve,
                              quality=args.quality_valve, limit=args.percentage_valve, start=args.start,
                              end=args.end, seq_size=args.seq_size)
    else:
        filtered1, filtered2 = filter_pe(fq1=f'"{args.fastq1}"', fq2=f'"{args.fastq2}"',
                                         o1=f'"{args.cleanq1}"', o2=f'"{args.cleanq2}"',
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
@arg_prop(dest='topology', choices=['linear', 'circular'], help=argparse.SUPPRESS, default='linear')
def annotate(args):
    # TODO:To fill the blanks of annotate method
    pass


@parse_func(func_help='visualization of GenBank file')
def visualize(args):
    # TODO:To fill the blanks of visualize method
    pass


@parse_func(func_help='run all the methods',
            parents=[universal_parser, assembly_parser, filter_parser, fastq_parser,
                     search_parser, saa_parser, annotation_parser])
@arg_prop(dest='disable_filter', help='filter will be not enabled if this switched on')
@arg_prop(dest='topology', choices=['linear', 'circular'], help=argparse.SUPPRESS, default='linear')
@arg_prop(dest='seq_size', help='how many sequences will be filtered out.', arg_type=int)
def all(args):

    # Go filtering
    if not args.disable_filter:
        args.fastq1, args.fastq2 = filter(args=args)

    contigs_file = assemble(args)

    # TODO:Finish findmitoscaf methods.


# This is ah, a somehow not ideal method in the whole MitoFlex coding,
# but it's done with the purpose of making the generated result
# cool and good.
# This acts as a cleaner, it cleans the generated empty folders (
# validated but not used) and some temporal files, it is executed no
# matter what the method is called, and even no matter if any the
# exposed function is called or NOT.
def cleanup(args):
    if args is None:
        return
    print(args)


# Entry starts at here
if __name__ == '__main__':

    parser = freeze_arguments('MitoFlex', desc)
    final_args = parse_then_call(parser)
    cleanup(final_args)
