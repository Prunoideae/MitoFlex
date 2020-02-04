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
    from utility.parser import parse_func, freeze_arguments, arg_prop, parse_then_call
    from utility.profiler import profiling
    from utility import logger
    # We are using this for making the main file clean, so wildcard is
    # not a problem here.
    from arguments import *  # pylint: disable=unused-wildcard-import

except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoFlex installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

# Debug only
logger.set_level(0)

# Constants
VERSION = '0.0.5'

# Command processing
desc = f"""
Description

    MitoFlex - A rewritten toolkit of its ancestor MitoZ for faster and better 
    mitochondrial assembly, annotation and visualization, for expandability and more.

Version
    {VERSION}

Citation
    MitoFlex
    
"""

@parse_func(func_help='filter out unqualified reads from fastq',
            parents=[universal_parser, fastq_parser, filter_parser])
def filter(args):

    if not path.isabs(args.cleanq1):
        args.cleanq1 = path.abspath(path.join(args.clean_dir, args.cleanq1))

    if args.cleanq2 is not None and not path.isabs(args.cleanq2):
        args.cleanq2 = path.abspath(path.join(args.clean_dir, args.cleanq2))

    filtered1 = filtered2 = None

    from filter.filter import filter_pe, filter_se

    if args.fastq2 is None:
        filtered1 = filter_se(fqiabs=f'"{args.fastq1}"', fqoabs=f'"{args.cleanq1}"', Ns=args.Ns_valve,
                              quality=args.quality_valve, limit=args.percentage_valve, start=args.start,
                              end=args.end)
    else:
        filtered1, filtered2 = filter_pe(fq1=args.fastq1, fq2=args.fastq2,
                                         o1=args.cleanq1, o2=args.cleanq2,
                                         dedup=args.deduplication,
                                         start=args.start, end=args.end,
                                         n=args.Ns_valve, q=args.quality_valve, l=args.percentage_valve)

    # Further processing for calling directly
    if args.__calling == 'filter':
        os.rename(filtered1, path.join(
            args.result_dir, path.basename(filtered1)))
        os.rename(filtered2, path.join(
            args.result_dir, path.basename(filtered2)))
    return filtered1, filtered2

@parse_func(func_help='assemble from input fastq reads, output contigs',
            parents=[universal_parser, fastq_parser, assembly_parser])
def assemble(args):

    from assemble.assemble import assemble as _assemble

    assembled_contigs = _assemble(fastq1=args.fastq1, fastq2=args.fastq2, base_dir=args.assemble_dir,
                                  work_prefix=args.workname, uselist=args.use_list,
                                  kmin=args.kmer_min, kmax=args.kmer_max, kstep=args.kmer_step, klist=args.kmer_list,
                                  no_mercy=args.no_mercy, disable_acc=args.disable_acc, prune_level=args.prune_level,
                                  prune_depth=args.prune_depth, keep_temp=not args.clean_temp, threads=args.threads)

    # Further processing for calling directly
    if args.__calling == 'assemble':
        os.rename(assembled_contigs, path.join(
            args.result_dir, path.basename(assembled_contigs)))

    return assembled_contigs


@parse_func(func_help='search for the most possible mitochondrial sequences from assembled data',
            parents=[universal_parser, fasta_parser, search_parser, saa_parser])
def findmitoscaf(args):

    from findmitoscaf.findmitoscaf import findmitoscaf as _findmitoscaf
    picked_fa = _findmitoscaf(
        thread_number=args.threads, clade=args.clade, relaxing=args.taxa_tolerance, gene_code=args.genetic_code,
        multi=args.min_abundance, taxa=args.required_taxa if not args.disable_taxa else None,
        prefix=args.workname, basedir=args.findmitoscaf_dir, contigs_file=args.fastafile, cover_valve=1)

    # Further processing for calling directly
    if args.__calling == 'findmitoscaf':
        os.rename(picked_fa, path.join(
            args.result_dir, path.basename(picked_fa)))
    return picked_fa


@parse_func(func_help='annotate PCGs, tRNA and rRNA genes',
            parents=[universal_parser, fasta_parser, annotation_parser, saa_parser, search_parser])
def annotate(args):
    from annotation.annotation import annotate as _annotate
    annotate_json, fa_file, rna_file = _annotate(basedir=args.annotation_dir, prefix=args.workname,
                                                 ident=30, fastafile=args.fastafile, genetic_code=args.genetic_code,
                                                 clade=args.clade, taxa=args.required_taxa, thread_number=args.threads,
                                                 wildcard_profile=args.wider_taxa)
    # Further processing for calling directly
    if args.__calling == 'annotate':
        import json
        with open(annotate_json, 'r') as f:
            pos_dict = json.load(f)
        with open(path.join(args.result_dir, 'annotated.txt'), 'w') as f:
            pcgs = {key: value for key, value in pos_dict.items()
                    if value[2] == 0}
            trna = {key: value for key, value in pos_dict.items()
                    if value[2] == 1}
            rrna = {key: value for key, value in pos_dict.items()
                    if value[2] == 2}

            print('PCGs found :')
            for key, value in pcgs.items():
                print(key, ':', value[0], '-', value[1], 'from', value[3])
            print('\ntRNAs found :')
            for key, value in trna.items():
                print(key, ':', value[0], '-', value[1], 'from', value[3])
            print('\nrRNAs found :')
            for key, value in rrna.items():
                print(key, ':', value[0], '-', value[1], 'from', value[3])
        os.rename(fa_file, path.join(args.result_dir, path.basename(fa_file)))
        os.rename(rna_file, path.join(
            args.result_dir, path.basename(rna_file)))

    return annotate_json


@parse_func(func_help='visualization of GenBank file')
def visualize(args):
    if not hasattr(args, 'use_json'):
        pass

    # TODO To fill the blanks of visualize method
    pass


@parse_func(func_help='run all the methods',
            parents=[universal_parser, assembly_parser, filter_parser, fastq_parser,
                     search_parser, saa_parser, annotation_parser])
@arg_prop(dest='disable_filter', help='filter will be not enabled if this switched on')
def all(args):

    # Go filtering
    if not args.disable_filter:
        #
        # Why I'm NOT using .gz ext here even I have implemented this:
        # 1. flate2 is slow, it takes much compressing data.
        # 2. plug in a SSD is much more easier than adding a CPU.
        #
        # You can still set this to xx.gz then it will surely make a
        # gzip for you, but this will have a great impact on the App's
        # running time, and it's strongly not recommended to do this.
        #
        args.cleanq1 = 'clean.1.fq'
        args.cleanq2 = 'clean.2.fq'
        args.fastq1, args.fastq2 = filter(args=args)

    args.fastafile = assemble(args)
    args.fastafile = findmitoscaf(args)

    if not args.disable_annotation:
        args.pos_json = annotate(args)
        args.use_json = True

        # Visualization is of no way if not annotated.
        # fastafile = findmitoscafed file
        # fastq1, fastq2 = filtered fastq file
        visualize(args)


def pre(args):
    # Initialize the logger.
    if hasattr(args, 'work_dir') and hasattr(args, 'workname'):
        logger.init(path.join(args.work_dir, f'{args.workname}.log'))
    else:
        logger.init(path.join(os.getcwd(), 'summary.log'))

    if hasattr(args, 'workname'):
        logger.log(2, f'MitoFlex {VERSION}, run {args.workname}')

    arg_dict = vars(args)
    logger.log(1, f'Arguments after parsed : ')
    logger.log(1, f'{[f"{key}={value}" for key, value in arg_dict.items()]}')

    if hasattr(args, 'disable_filter') and args.disable_filter:
        logger.log(3, 'Filtering is not enabled.')

    if hasattr(args, 'disable_annotation') and args.disable_annotation:
        logger.log(3, 'Annotation is not enabled.')


def post(args):
    if args is None:
        return
    logger.finalize()


# Entry starts at here
if __name__ == '__main__':
    parser = freeze_arguments('MitoFlex', desc)
    parse_then_call(parser, pre=pre, post=post)
