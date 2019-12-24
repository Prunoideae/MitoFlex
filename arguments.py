"""
arguments.py
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

import sys
import os

try:
    from utility.parser import register_group
except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoX installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

'''

'''
# Universal group


def universal_regulator(args):
    tr = os.cpu_count()
    tr = tr if tr is not None else 8
    if args.threads <= 0 or args.threads > tr:
        print(
            f"Specified thread number not in range, using {tr} threads instead.")
        args.threads = tr

    args.work_dir = os.path.abspath(os.path.join(args.basedir, args.workname))
    args.result_dir = os.path.abspath(os.path.join(
        args.work_dir, args.workname + '.result'))
    args.temp_dir = os.path.abspath(os.path.join(
        args.work_dir, args.workname + '.temp'))
    if args.profiling:
        args.profile_dir = os.path.abspath(os.path.join(
            args.temp_dir, 'performance'))

    # Validates the folders
    try:
        os.makedirs(args.work_dir, exist_ok=True)
        os.makedirs(args.result_dir, exist_ok=True)
        os.makedirs(args.temp_dir, exist_ok=True)
        if args.profiling:
            os.makedirs(args.profile_dir, exist_ok=True)
    except Exception as identifier:
        print(
            'Error occured when validating the directories, please check your permissions or things could be related.')
        return False
    return True


universal_parser, universal_group = register_group('Universal arguments', [
    {
        'name': 'workname',
        'required': True,
        'help': 'work output name.'
    },
    {
        'name': 'threads',
        'default': 8,
        'help': 'thread numbers.'
    },
    {
        'name': 'clean-temp',
        'default': False,
        'help': 'remove temporal files and folders after work done.'
    },
    {
        'name': 'basedir',
        'default': os.getcwd(),
        'help': 'working folder will be generated in which directory.'
    },
    {
        'name': 'profiling',
        'default': False,
        'help': 'a performance profiling will be generated if this is on.'
    }], func=universal_regulator
)

# Fastq group


def fastq_regulator(args):
    valid = True

    if args.fastq1 is None and args.fastq2 is not None:
        args.fastq1, args.fastq2 = args.fastq2, args.fastq1

    if args.fastq1 is None and args.fastq2 is None:
        print("Both fastq file inputs are missing. Exiting the run.")
        return False

    # Explicitly using cwd path here.
    args.fastq1 = os.path.abspath(args.fastq1)
    if args.fastq2 is not None:
        args.fastq2 = os.path.abspath(args.fastq2)

    if not (os.path.isfile(args.fastq1) and (args.fastq2 is None or os.path.isfile(args.fastq2))):
        valid = False
        print("Input FASTQ file is not valid.")

    if args.fastq_read_length <= 0:
        valid = False
        print("Specified fastq read length is not valid.")

    if hasattr(args, 'use_list'):
        min_kmer = int(args.kmer_list.split(
            ',')[0]) if args.use_list else args.kmer_min
        if min_kmer >= args.fastq_read_length:
            valid = False
            print("Specified fastq read length lower than the mininum kmer.")

    return valid


fastq_parser, fastq_group = register_group('Fastq arguments', [
    {
        'name': 'fastq1',
        'meta': 'file',
        'help': 'fastq file 1 input.'
    },
    {
        'name': 'fastq2',
        'meta': 'file',
        'help': 'fastq file 2 input.'
    },
    {
        'name': 'fastq-alter-format',
        'default': False,
        'help': 'using FASTQ-like (Q+64) format or Q+33 file format.'
    },
    {
        'name': 'fastq-read-length',
        'default': 150,
        'help': '''read length of fastq reads, used in MEGAHIT and bwa,
        it must be at least larger than the k-min if assemble specified.'''
    }
], func=fastq_regulator)

# Fasta group


def fasta_regulator(args):
    args.fastafile = os.path.abspath(args.fastafile)
    if not os.path.isfile(args.fastafile):
        print("Input FASTA file not valid.")
        return False
    return True


fasta_parser, fasta_group = register_group('Fasta arugments', [
    {
        'name': 'fastafile',
        'meta': 'file',
        'help': 'fasta file'
    }
], func=fasta_regulator)

# Filter group


def filter_regulator(args):

    if hasattr(args, 'disable_filter') and args.disable_filter:
        return True

    valid = True

    try:
        args.start, args.end, *_ = [int(x) if int(x) > -1 else None
                                    for x in args.keep_region.split(',')]
    except Exception as identifier:
        print('Input range is not valid.')
        valid = False

    if args.adapter1 is not None:
        args.adapter1 = os.path.abspath(args.adapter1)
    if args.adapter2 is not None:
        args.adapter2 = os.path.abspath(args.adapter2)

    if hasattr(args, 'temp_dir'):
        args.clean_dir = os.path.join(args.temp_dir, 'cleandata')
    else:
        args.clean_dir = os.getcwd()

    try:
        os.makedirs(args.clean_dir, exist_ok=True)
    except Exception as identifier:
        valid = False
        print(
            'Error occured when validating the directories, please check your permissions or things could be related.')

    if not ((args.adapter1 is None or os.path.isfile(args.adapter1) and (args.adapter2 is None or os.path.isfile(args.adapter2)))):
        print('Input adapter file is not valid.')
        valid = False

    if args.quality_valve <= 0 or args.quality_valve >= 255:
        print('Input quality limit is not valid.')
        valid = False

    if args.Ns_valve <= 0:
        print('Input N limit is not valid.')
        valid = False

    if args.percentage_valve <= 0 or args.percentage_valve >= 1:
        print('Input percentage limit is not valid.')
        valid = False

    return valid


filter_parser, filter_group = register_group('Filter argumetns', [
    {
        'name': 'cleanq1',
        'meta': 'file',
        'help': 'cleandata output file 1'
    },
    {
        'name': 'cleanq2',
        'meta': 'file',
        'help': 'cleandata output file 2'
    },
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
], func=filter_regulator)

# Assembly group


def assembly_regulator(args):
    valid = True

    if args.insert_size <= 0:
        valid = False
        print('Input insert size is not valid.')

    if hasattr(args, 'temp_dir'):
        args.assemble_dir = os.path.join(args.temp_dir, 'assemble')
    else:
        args.assemble_dir = os.getcwd()

    try:
        os.makedirs(args.assemble_dir, exist_ok=True)
    except Exception as identifier:
        valid = False
        print(
            'Error occured when validating the directories, please check your permissions or things could be related.')

    if args.use_list:
        args.kmer_list = [int(x) for x in args.kmer_list.split(',')]
        args.kmer_list.sort()
        if 0 in [x % 2 for x in args.kmer_list]:
            print('All kmer length must be odd.')
            valid = False
    else:
        if True in [
            args.kmer_min <= 0,
            args.kmer_max <= 0,
            args.kmer_step <= 0,
            args.kmer_max < args.kmer_min
        ]:
            print('Input kmer arguments have invalid values.')
            valid = False

        if args.kmer_min % 2 == 0 or (args.kmer_min + args.kmer_step) % 2 == 0:
            print('All kmer length must be odd.')
            valid = False

    if args.prune_depth < 0:
        print('Prune depth lower than 0.')

    return valid


assembly_parser, assembly_group = register_group('Assembly arguments', [
    {
        'name': 'insert-size',
        'default': 150,
        'help': 'insert size of input fastq files.'
    },
    {
        'name': 'use-list',
        'default': False,
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
        'help': 'force disable the hardware accerlation if there\'s some problem in using it, only use as a last resort.'
    }
], func=assembly_regulator)

# Search mitochondrial gene group


def search_regulator(args):
    valid = True

    if args.min_abundance <= 0:
        print("Input minimum abundance is not valid.")
        valid = False

    from ete3 import NCBITaxa
    ncbi = NCBITaxa()

    if hasattr(args, 'temp_dir'):
        args.assemble_dir = os.path.join(args.temp_dir, 'assemble')
    else:
        args.assemble_dir = os.getcwd()

    try:
        os.makedirs(args.assemble_dir, exist_ok=True)
    except Exception as identifier:
        valid = False
        print(
            'Error occured when validating the directories, please check your permissions or things could be related.')

    if args.required_taxa not in ncbi.get_name_translator([args.required_taxa]):
        print("Specified taxanomy name not in NCBI taxanomy database.")
        return False
    args.taxa_ids = ncbi.get_name_translator([args.required_taxa])[
        args.required_taxa]

    return valid


search_parser, search_group = register_group('Search mitochondrial sequences arguments', [
    {
        'name': 'filter-taxa',
        'default': False,
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
        'help': 'taxa sequences other than this will not be selected.'
    },
    {
        'name': 'taxa-tolerance',
        'default': 0,
        'choices': list(range(7)),  # 0-6
        'help': 'a tolerance for non-target sequences.'
    }
], func=search_regulator)

# Search and annotation mitochondrial gene group


def saa_regulator(args):
    gene_code = {
        # Standard Code 1, raise a warning because it's not commonly used.
        '#1': 'Standard code',

        # Explicit mitochondrial code.
        '2': 'Vertebrate Mitochondrial code',
        '3': 'Yeast Mitochondrial code',
        '4': 'Mold, Protozoan and Coelenterate Mitochondiral code and the Mycoplasma/Spiroplasma code',
        '5': 'Invertebrate Mitochondrial code',
        '9': 'Echinoderm and Flatworm Mitochondrial code',
        '13': 'Ascidian Mitochondrial code',
        '14': 'Alternative Flatworm Mitochondrial code',
        '16': 'Chlorophycean Mitochondrial code',
        '21': 'Trematode Mitochondrial code',
        '22': 'Scenedesmus obliquus Mitochondrial code',
        '23': 'Thraustochytrium Mitochondrial code',
        '24': 'Pterobranchia Mitochondrial code',
        '33': 'Cephalodiscidae Mitochondrial UAA-Tyr code',

        # Prokaryotic or Plastid code, THEORETICALLY correct, but untested.
        '_11': 'Baterial, Archaeal and Plant Plastid code',
        '_25': 'Candidate Division SR1 and Gracilibacteria code',

        # Nuclear code, untested and with high risk of the whole toolkit.
        '-6': 'Ciliate, Dasycladacean and Hexamita Nuclear code',
        '-10': 'Euplotid Nuclear code',
        '-12': 'Alternative Yeast Nuclear code',
        '-26': 'Pachysolen tannophilus Nuclear Code',
        '-27': 'Karyorelict Nuclear code',
        '-28': 'Condylostoma Nuclear code',
        '-29': 'Mesodinium Nuclear code',
        '-30': 'Peritrich Nuclear code',
        '-31': 'Blastocrithidia Nuclear code'
    }

    code = str(args.genetic_code)

    if code in gene_code:
        print(f'Using genetic code {args.genetic_code} : {gene_code[code]}')
    elif '_' + code in gene_code:
        print(
            f'Using genetic code {args.genetic_code} : {gene_code["_" + code]}.\nThese codes are thoretically suitable with current workflow, but it\'s NOT tested.')
    elif '-' + code in gene_code:
        print(
            f'Using genetic code {args.genetic_code} : {gene_code["-" + code]}!\nThey are obviously out of MitoX\'s range, so further calls needs your validation.')
        # Only args with recessive option y could bypass this, thus every dangerous method should be handled by hand unless a config say so.
        if not hasattr(args, 'y') or not args.y:
            answer = input("Continue? Y/[N] : ")
            while answer not in ['Y', 'N', '', 'y', 'n']:
                answer = input("Continue? Y/[N] : ")
            if answer.upper() is not 'Y':
                sys.exit('Exited.')
    elif code == '1':
        print('You are using Standard Code 1! Please make sure you really have the need to do this, because the standard code is NOT meant to be applied in most cases!')
        if not hasattr(args, 'y') or not args.y:
            answer = input("Continue? Y/[N] : ")
            while answer not in ['Y', 'N', '', 'y', 'n']:
                answer = input("Continue? Y/[N] : ")
            if answer.upper() is not 'Y':
                sys.exit('Exited.')
    else:
        print('Genetic code not found in the NCBI table!')
        return False
    return True


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
], func=saa_regulator)

# Annotation group
annotation_parser, annotation_group = register_group('Annotation arguments', [
    {
        'name': 'disable-annotation',
        'default': True,
        'help': 'switching this on will disable annotation.'
    },
    {
        'name': 'species-name',
        'default': 'Test sp.',
        'help': 'species name to use in genbank file.'
    }
])
