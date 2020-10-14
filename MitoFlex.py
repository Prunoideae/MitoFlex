#!/usr/bin/env python3

"""
MitoFlex
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

from os import path
import os
import sys
import traceback
import time

if sys.version_info[0] < 3:
    sys.exit('Python 3 must be installed in current environment! Please check if any of your environment setup (like conda environment) is deactivated or wrong!')

try:
    from utility.parser import parse_func, freeze_arguments, arg_prop, parse_then_call
    from utility import logger
    from utility.helper import shell_call
    from arguments import *  # pylint: disable=unused-wildcard-import
    import configurations
except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoFlex installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

# Constants
VERSION = '0.2.9'

# Static variables
start_time = time.time()

# Command processing
desc = f"""
Description
    An almost all-in-one pipeline for Mitogenome analysis from de novo NGS data
    aiming at result quality, speed and flexibility.

Version
    {VERSION}

Citation
    MitoFlex: an efficient, high-performance toolkit for mitogenome assembly, 
    annotation, and visualization. (Not yet published)

"""

@parse_func(func_help='filter out unqualified reads from fastq',
            parents=[universal_parser, fastq_parser, filter_parser])
def filter(args):

    if not hasattr(args, 'disable_filter'):
        args.disable_filter = False

    if not args.cleanq1:
        args.cleanq1 = 'clean.1.fq'
    if args.fastq2 and not args.cleanq2:
        args.cleanq2 = 'clean.2.fq'

    args.cleanq1 = path.abspath(path.join(args.clean_dir, args.cleanq1))
    if hasattr(args, 'fastq2') and args.fastq2:
        args.cleanq2 = path.abspath(path.join(args.clean_dir, args.cleanq2))

    filtered1 = filtered2 = None

    from filter.filter import filter_pe, filter_se

    if args.fastq2 is None:
        filtered1 = filter_se(fqiabs=args.fastq1, fqoabs=args.cleanq1, Ns=args.Ns_valve,
                              quality=args.quality_valve, limit=args.percentage_valve, start=args.start,
                              end=args.end, trim=args.trimming, trunc=args.disable_filter)
    else:
        filtered1, filtered2 = filter_pe(fq1=args.fastq1, fq2=args.fastq2,
                                         o1=args.cleanq1, o2=args.cleanq2,
                                         dedup=args.deduplication,
                                         start=args.start, end=args.end,
                                         n=args.Ns_valve, q=args.quality_valve, l=args.percentage_valve, trim=args.trimming,
                                         trunc=args.disable_filter)

    # Further processing for calling directly
    if args.__calling == 'filter':
        os.rename(filtered1, path.join(
            args.result_dir, path.basename(filtered1)))
        if filtered2:
            os.rename(filtered2, path.join(
                args.result_dir, path.basename(filtered2)))
    return filtered1, filtered2


@parse_func(func_help='assemble from input fastq reads, output contigs',
            parents=[universal_parser, fastq_parser, assembly_parser])
def assemble(args):

    from assemble.assemble import assemble as _assemble

    assembled_contigs = _assemble(fastq1=args.fastq1, fastq2=args.fastq2, base_dir=args.assemble_dir,
                                  work_prefix=args.workname, disable_local=args.disable_local,
                                  kmer_list=args.kmer_list, depth_list=args.depth_list,
                                  prune_level=args.prune_level, prune_depth=args.prune_depth,
                                  keep_temp=args.keep_temp, threads=args.threads,
                                  insert_size=args.insert_size, no_scaf=args.disable_scaffolding)

    # Further processing for calling directly
    if args.__calling == 'assemble':
        os.rename(assembled_contigs, path.join(
            args.result_dir, path.basename(assembled_contigs)))

    return assembled_contigs


@parse_func(func_help='search for the most possible mitochondrial sequences from assembled data',
            parents=[universal_parser, fasta_parser, search_parser, saa_parser, fastq_parser])
@arg_prop(dest='from_megahit', help='on if the result is from megahit, so remapping will be skipped.', default=False)
def findmitoscaf(args):

    if args.__calling == 'findmitoscaf':

        if not args.from_megahit:
            logger.log(2, 'Remapping reads to contigs since contigs are not assembled from pipeline.')
            fastfilter_bin = path.abspath(path.join(path.dirname(__file__), 'assemble', 'fastfilter'))
            filtered_fasta = path.join(args.findmitoscaf_dir, f'{args.workname}.filtered.fa')
            shell_call(fastfilter_bin, i=args.fastafile, o=filtered_fasta,
                       l=f"{configurations.assemble.min_length},{configurations.assemble.max_length}",
                       d=0)
            fq1, fq2 = args.fastq1, args.fastq2
            if not (fq1 or fq2):
                raise RuntimeError("At least one fastq file should be specified!")
            if not fq1:
                fq1, fq2 = fq2, fq1
            # Remapping to calculate average depth.
            from findmitoscaf.findmitoscaf import remap_sequence
            args.fastafile = remap_sequence(args.workname, args.findmitoscaf_dir, filtered_fasta, args.fastq1, args.fastq2, args.threads)
        else:
            logger.log(2, "Remapping skipped since from-megahit is specified, no tagging needed.")

    from findmitoscaf.findmitoscaf import findmitoscaf as _findmitoscaf
    picked_fa = _findmitoscaf(
        thread_number=args.threads, clade=args.clade, relaxing=args.taxa_tolerance, gene_code=args.genetic_code,
        multi=args.min_abundance, taxa=args.required_taxa if not args.disable_taxa else None,
        prefix=args.workname, basedir=args.findmitoscaf_dir, contigs_file=args.fastafile,
        merge_method=args.merge_method, merge_overlapping=args.merge_overlap, merge_search=args.merge_start)

    # Further processing for calling directly
    if args.__calling == 'findmitoscaf':
        os.rename(picked_fa, path.join(
            args.result_dir, path.basename(picked_fa)))
    return picked_fa


@parse_func(func_help='annotate PCGs, tRNA and rRNA genes',
            parents=[universal_parser, fasta_parser, annotation_parser, saa_parser, search_parser])
def annotate(args):

    from annotation.annotation import annotate as _annotate, fix_circular

    # Check assemble file, if only one sequence and itself is circular, the genome is then circular.
    circular = False
    if configurations.annotation.trim_circular:
        circular = fix_circular(fa_file=args.fastafile)

    # Annotate the file
    annotate_json, fa_file, rna_file = _annotate(basedir=args.annotation_dir, prefix=args.workname,
                                                 ident=30, fastafile=args.fastafile, genetic_code=args.genetic_code,
                                                 clade=args.clade, thread_number=args.threads,
                                                 wildcard_profile=args.wider_taxa, trna_overlapping=30,
                                                 hmmer_search=args.use_hmmer, score=args.hmmer_score, e_value=args.hmmer_e)

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
        if circular:
            print('The final mitogenome is circular and trimmed.')
        os.rename(fa_file, path.join(args.result_dir, path.basename(fa_file)))
        os.rename(rna_file, path.join(
            args.result_dir, path.basename(rna_file)))

    return annotate_json, circular, fa_file, rna_file


@parse_func(func_help='visualization of sequences',
            parents=[universal_parser, fasta_parser, fastq_parser])
@arg_prop(dest='pos_json', help='specify the json file for marking genes')
@arg_prop(dest='circular', help='draw the genome like a circle or have some break', default=False)
def visualize(args):

    basedir = args.temp_dir if args.__calling == 'visualize' else path.join(
        args.temp_dir, 'visualize')
    try:
        os.makedirs(basedir, exist_ok=True)
    except Exception:
        raise RuntimeError("Cannot validate folder for visualization!")

    from visualize.visualize import visualize as _visualize
    circos_png, circos_svg = _visualize(fasta_file=args.fastafile, fastq1=args.fastq1, fastq2=args.fastq2,
                                        pos_json=args.pos_json, prefix=args.workname, basedir=basedir,
                                        threads=args.threads, circular=args.circular)

    # Further processing for calling directly
    if args.__calling == 'visualize':
        os.rename(circos_png, path.join(
            args.result_dir, path.basename(circos_png)))
        os.rename(circos_svg, path.join(
            args.result_dir, path.basename(circos_svg)))
        return None, None

    return circos_png, circos_svg


@parse_func(func_help='run all the methods',
            parents=[universal_parser, assembly_parser, filter_parser, fastq_parser,
                     search_parser, saa_parser, annotation_parser])
@arg_prop(dest='disable_filter', help='filter will be not enabled if this switched on', default=False)
@arg_prop(dest='disable_visualization', help='visualization will be not enabled if this switched on', default=False)
def all(args):

    # Go filtering
    #
    # Why I'm NOT using .gz ext here even I have implemented this:
    # 1. flate2 is slow, it takes much compressing data if single-threaded.
    # 2. plug in a SSD is much more easier than adding a CPU.
    # 3. Some method uses only plain text data, so you need an extra (de)compression
    #    but it means nothing in the process.
    # 4. Some further codes may only accept plain-text input, and I'm not adding
    #    support of gzip to it.

    args.cleanq1 = 'clean.1.fq'
    args.cleanq2 = 'clean.2.fq'
    if configurations.filter_rawdata.compress_output_in_all:
        args.cleanq1 += '.gz'
        args.cleanq2 += '.gz'

    args.fastq1, args.fastq2 = filter(args)

    args.fastafile = assemble(args)
    args.fastafile = findmitoscaf(args)

    if not args.disable_annotation:
        (args.pos_json, args.circular,
         args.annotated_cds, args.annotated_rna) = annotate(args)

        # Visualization is of no way if not annotated.
        args.circos_png, args.circos_svg = visualize(
            args) if not args.disable_visualization else (None, None)

    # Add command check if there's something further
    # If you wrapped the 'all' module in other task or workflow
    # the results will be retained since we don't know what you
    # want.
    if args.__calling == 'all':
        def move_to_result(*files):
            for file in files:
                if path.isfile(str(file)):
                    os.rename(file, path.join(
                        args.result_dir, path.basename(file)))
        # Iteratively collects all the results generated in the whole process
        move_to_result(args.circos_png, args.circos_svg,
                       args.pos_json, args.fastafile,
                       args.annotated_cds, args.annotated_rna)
        logger.log(2, f'Results dumped at {args.result_dir}')


@parse_func(func_help='load all modules provided by MitoFlex, use to test if some modules are not installed correctly.')
def load_modules(args):
    try:
        logger.log(2, 'Loading filter module.')
        from filter.filter import filter_pe, filter_se
        logger.log(2, 'Loading assemble module.')
        from assemble.assemble import assemble
        logger.log(2, 'Loading findmitoscaf module.')
        from findmitoscaf.findmitoscaf import findmitoscaf
        logger.log(2, 'Loading annotation module.')
        from annotation.annotation import annotate
        logger.log(2, 'Loading visualize module.')
        from visualize.visualize import visualize
    except Exception:
        logger.log(4, 'Cannot load module!')
    finally:
        logger.log(2, 'All modules are loaded correctly.')


# This is for initializing the framework right before the command executed,
# but after the arguments are processed. Pre will initialize something no
# matter what command is called. Not pretty.
def pre(args):

    # Initialize the logger.
    if hasattr(args, 'work_dir') and hasattr(args, 'workname'):
        logger.init(path.join(args.work_dir, f'{args.workname}.log'))
    else:
        logger.init(path.join(os.getcwd(), 'summary.log'))
    if hasattr(args, 'level'):
        logger.set_level(args.level)
    logger.log(
        2, f'MitoFlex {VERSION}, run {args.workname if hasattr(args, "workname") else "1"}')

    arg_dict = vars(args)
    logger.log(2, f'Arguments after parsed : ')
    logger.log(2, f'{[f"{key}={value}" for key, value in arg_dict.items()]}')

    if hasattr(args, 'disable_filter') and args.disable_filter:
        logger.log(3, 'Filtering is not enabled, files will only be truncated.')

    if hasattr(args, 'disable_annotation') and args.disable_annotation:
        logger.log(3, 'Annotation is not enabled.')

    def runtime_error_logger(exception_type, value, tb):
        if exception_type == RuntimeError:
            logger.log(4, value)
            logger.log(
                4,
                'A RuntimeError was occured! This is already considered in the code'
                ', but since it\'s thought to be errors in parts outside the MitoFlex can handle, it\'s'
                ' NOT a bug caused by MitoFlex itself. Please check the error message'
                ' and try to fix the possible cause of the crash, only as a last resort, send '
                'github a issue with a rerun with logger level set to 0.'
            )
            logger.finalize()
            sys.exit()
        else:
            if exception_type != KeyboardInterrupt:
                logger.log(
                    4,
                    "An unexpected error was happened in the MitoFlex, this could be an bug in the program,"
                    " so please report it if you see this message in log.")
                logger.log(
                    4, f"Error type : {exception_type.__name__}, value : {value}")
                logger.log(
                    4, f"Traceback :")
                logger.__log(traceback.format_tb(tb=tb))
                with open(path.join(path.dirname(logger.get_file()), 'traceback.txt'), 'w') as f:
                    traceback.print_tb(tb, file=f)

                logger.log(4, "Logging additional information")
                import psutil
                curp = psutil.Process()
                logger.log(4, curp.open_files())
                logger.log(4, curp.environ())
                logger.log(4, curp.memory_full_info())

            else:
                logger.log(2, "This run was terminated manually.")
            logger.finalize()
            sys.__excepthook__(exception_type, value, tb)

    sys.excepthook = runtime_error_logger


# This is for cleaning up temporal files generated by commands, and clean up
# the environment to make a proper end.
def post(args):

    if args is None:
        return
    if hasattr(args, 'keep_temp') and not args.keep_temp and args.__calling != 'filter' and hasattr(args, 'cleanq1'):
        # Not removing until here since cleanq1 and cleanq2 have many other usage other than assembling
        logger.log(1, 'Removing filtered data files.')
        os.remove(args.cleanq1)
        if args.fastq2 != None:
            os.remove(args.cleanq2)
    logger.log(2, f'All done! Time elapsed : {time.time()-start_time:.2f}s.')
    logger.finalize()


# Entry starts at here, things are mainly handled by the built-in CLI, please
# do not change anything here unless you really know what you are doing.
if __name__ == '__main__':
    parser = freeze_arguments('MitoFlex', desc)
    parse_then_call(parser, pre=pre, post=post)
