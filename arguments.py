"""
arguments.py
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

import sys
import os

try:
    from utility.parser import register_group
except ModuleNotFoundError as identifier:
    print(
        f'Module {identifier.name} not found! Please check your MitoFlex installation!')
    sys.exit()
except ImportError as identifier:
    print(
        f'Error occured when importing module {identifier.name}! Please check your system, python or package installation!')
    sys.exit()

profile_dir_tbn = os.path.join(
    os.path.dirname(__file__), 'profile', 'MT_database')
profile_dir_hmm = os.path.join(os.path.dirname(__file__), 'profile', 'CDS_HMM')


def safe_makedirs(path, exist_ok=False):
    if os.path.isdir(path) and not exist_ok:
        raise FileExistsError()
    os.makedirs(path, exist_ok=True)


# Universal group
def universal_regulator(args):

    args.work_dir = os.path.abspath(os.path.join(args.basedir, args.workname))

    args.result_dir = os.path.abspath(os.path.join(
        args.work_dir, args.workname + '.result'))
    args.temp_dir = os.path.abspath(os.path.join(
        args.work_dir, args.workname + '.temp'))
    # Validates the folders
    try:
        safe_makedirs(args.work_dir)
        safe_makedirs(args.result_dir)
        safe_makedirs(args.temp_dir)
    except FileExistsError:
        print('Diretory exist before validating, please check and remove it to prevent data loss.')
        return False
    except Exception:
        print(
            'Error occured when validating the directories, please check your permissions or things could be related.')
        return False

    tr = os.cpu_count()
    tr = tr if tr is not None else 8
    if args.threads <= 0 or args.threads > tr:
        print(
            f"Specified thread number not in range, using {tr} threads instead.")
        args.threads = tr

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
        'name': 'keep-temp',
        'default': False,
        'help': 'remove temporal files and folders after work done.'
    },
    {
        'name': 'basedir',
        'default': os.getcwd(),
        'help': 'working folder will be generated in which directory.'
    },
    {
        'name': 'level',
        'default': 2,
        'choices': [0, 1, 2, 3, 4],
        'help': 'how verbose will MitoFlex output logs, the lower the log outputs more, level 0 is mainly for code debugging and reporting.'
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

    if args.level == 0:
        print("Using logger level 0, this could lead to a very verbose logging.")

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
        'help': 'fasta file, automatically assembled in all method'
    }
], func=fasta_regulator)


# Filter group
def filter_regulator(args):

    if hasattr(args, 'disable_filter') and args.disable_filter:
        return True

    valid = True

    try:
        args.start, args.end, *_ = [int(x) if int(x) > 0 else 0
                                    for x in args.keep_region.split(',')]
    except Exception:
        print('Input range is not valid.')
        valid = False

    if hasattr(args, 'temp_dir'):
        args.clean_dir = os.path.join(args.temp_dir, 'cleandata')
    else:
        args.clean_dir = os.path.join(os.getcwd(), 'cleandata')

    try:
        safe_makedirs(args.clean_dir)
    except FileExistsError:
        valid = False
        print('Diretory exist before validating, please check and remove it to prevent data loss.')
    except Exception:
        valid = False
        print('Error occured when validating the directories, please check your permissions or things could be related.')

    if args.quality_valve <= 0 or args.quality_valve >= 255:
        print('Input quality limit is not valid.')
        valid = False

    if args.Ns_valve <= 0:
        print('Input N limit is not valid.')
        valid = False

    if args.percentage_valve <= 0 or args.percentage_valve >= 1:
        print('Input percentage limit is not valid.')
        valid = False
    if args.trimming < 0:
        print('Trimming value is not valid.')
        valid = False
    args.trimming = int(float(args.trimming) * (10**9))
    return valid


filter_parser, filter_group = register_group('Filter arguments', [
    {
        'name': 'cleanq1',
        'meta': 'file',
        'help': 'cleandata output file 1, automatically set in all method'
    },
    {
        'name': 'cleanq2',
        'meta': 'file',
        'help': 'cleandata output file 2, automatically set in all method'
    },
    {
        'name': 'deduplication',
        'default': False,
        'help': 'fitler duplication caused by adapter ligation if switched on.'
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
        'default': '0,0',
        'meta': 'beg,end',
        'help': 'only the BEG and the END will be read, leave blank for full length.'
    },
    {
        'name': 'trimming',
        'default': 5,
        'help': 'only filter out x Gbps of the given dataset.'
    }
], func=filter_regulator)


# Assembly group
def assembly_regulator(args):
    valid = True

    if hasattr(args, 'temp_dir'):
        args.assemble_dir = os.path.join(args.temp_dir, 'assemble')
    else:
        args.assemble_dir = os.path.join(os.getcwd(), 'assemble')

    try:
        safe_makedirs(args.assemble_dir)
    except FileExistsError:
        print('Diretory exist before validating, please check and remove it to prevent data loss.')
        valid = False
    except Exception:
        valid = False
        print('Error occured when validating the directories, please check your permissions or things could be related.')

    if args.use_list:
        args.kmer_list = [int(x) for x in args.kmer_list.split(',')]
        args.kmer_list.sort()
        if 0 in [x % 2 for x in args.kmer_list]:
            print('All kmer length must be odd.')
            valid = False
        args.kmer_list = ','.join([str(x) for x in args.kmer_list])
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
        valid = False

    try:
        args.additional_kmers = [
            int(x)
            for x in args.additional_kmers.split(',')
            if x
        ]
    except Exception:
        print('Error parsing final kmers.')
        valid = False

    return valid


assembly_parser, assembly_group = register_group('Assembly arguments', [
    {
        'name': 'use-list',
        'default': False,
        'help': 'k-mer list arguments will be used if this switched on.'
    },
    {
        'name': 'disable-local',
        'default': False,
        'help': 'disable the local assembly of MEGAHIT, will have a performance imporve but may lowers the result quality.'
    },
    {
        'name': 'kmer-list',
        'default': '31,39,59,79,99,119,141',
        'help': 'list of kmer to use in sDBG building, all length must be odd.'
    },
    {
        'name': 'kmer-min',
        'default': 31,
        'help': 'the minimum length kmer used in MEGAHIT\'s multiple kmer assemble strategy.'
    },
    {
        'name': 'kmer-max',
        'default': 141,
        'help': 'the maximum length kmer used in MEGAHIT\'s multiple kmer assemble strategy.'
    },
    {
        'name': 'kmer-step',
        'default': 22,
        'help': 'increment of kmer size of each iteration (<= 28), must be even number.'
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
        'name': 'additional-kmers',
        'default': '',
        'help': 'input a list of kmers seperated by comma, to specify what kmer results will be added in the end.'
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
        args.findmitoscaf_dir = os.path.join(args.temp_dir, 'findmitoscaf')
    else:
        args.findmitoscaf_dir = os.path.join(os.getcwd(), 'findmitoscaf')
    try:
        safe_makedirs(args.findmitoscaf_dir)
    except FileExistsError:
        valid = False
        print('Diretory exist before validating, please check and remove it to prevent data loss.')
    except Exception:
        valid = False
        print('Error occured when validating the directories, please check your permissions or things could be related.')

    if args.required_taxa not in ncbi.get_name_translator([args.required_taxa]):
        print("Specified taxanomy name not in NCBI taxanomy database.")
        return False
    args.taxa_ids = ncbi.get_name_translator([args.required_taxa])[
        args.required_taxa]

    return valid


search_parser, search_group = register_group('Search mitochondrial sequences arguments', [
    {
        'name': 'disable-taxa',
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
        'default': 'Arthropoda',
        'help': 'taxa sequences other than this taxanomy class will not be selected.'
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

        print(
            f'Using genetic code {args.genetic_code} : {gene_code[code]}')
    elif '_' + code in gene_code:

        print(
            f'Using genetic code {args.genetic_code} : {gene_code["_" + code]}.\nThese codes are thoretically suitable with current workflow, but it\'s NOT tested.')
    elif '-' + code in gene_code:
        print(
            f'Using genetic code {args.genetic_code} : {gene_code["-" + code]}!\nThey are obviously out of MitoFlex\'s range, so further calls needs your validation.')
        # Only args with recessive option y could bypass this, thus every dangerous method should be handled by hand unless a config say so.
        if not hasattr(args, 'y') or not args.y:
            answer = input("Continue? Y/[N] : ")
            while answer not in ['Y', 'N', '', 'y', 'n']:
                answer = input("Continue? Y/[N] : ")
            if answer.upper() is not 'Y':
                sys.exit('Exited.')

            print(
                f'Using genetic code {code} : {gene_code["-" + code]}, assembling quality may not be guaranteed.')
    elif code == '1':
        print('You are using Standard Code 1! Please make sure you really have the need to do this, because the standard code is NOT meant to be applied in most cases!')
        if not hasattr(args, 'y') or not args.y:
            answer = input("Continue? Y/[N] : ")
            while answer not in ['Y', 'N', '', 'y', 'n']:
                answer = input("Continue? Y/[N] : ")
            if answer.upper() is not 'Y':
                sys.exit('Exited.')

            print(f'Using genetic code 1, assembling quality can\'t be guaranteed!')
    else:
        print('Genetic code not found in the NCBI table!')
        return False
    return True


profiles_hmm = [os.path.splitext(profile_name)[0]
                for profile_name in os.listdir(profile_dir_hmm)]
saa_parser, saa_group = register_group('Search and annotate arguments', [
    {
        'name': 'genetic-code',
        'default': 5,
        'help': 'genetic code table to be used in the run.'
    },
    {
        'name': 'clade',
        'default': 'Arthropoda',
        'choices': [os.path.splitext(profile_name)[0]
                    for profile_name in os.listdir(profile_dir_tbn)
                    if os.path.splitext(profile_name)[0] in profiles_hmm],
        'help': 'which clade\'s nhmmer profile and cds database will be used in the run.'
    },
    {
        'name': 'max-contig-length',
        'default': 20000,
        'help': 'the maximum length of a reasonable sequence. Use to filter out strangly connected scaffolds caused by containmination or other reasons before picking.'
    },
    {
        'name': 'wider-taxa',
        'default': False,
        'help': 'to use the general protein profile for annotation, only used if no reports and cannot add addtional profile.'
    },



], func=saa_regulator)


def annotation_regulator(args):
    if hasattr(args, 'temp_dir'):
        args.annotation_dir = os.path.join(args.temp_dir, 'annotation')
    else:
        args.annotation_dir = os.path.join(os.getcwd(), 'annotation')

    try:
        safe_makedirs(args.annotation_dir)
    except FileExistsError:
        print('Diretory exist before validating, please check and remove it to prevent data loss.')
        return False
    except Exception:
        print('Cannot make annotation temp folder!')
        return False

    valid = True

    if args.max_contig_length < 0:
        print('Please input valid contig length.')
        valid = False
    if args.hmmer_score < 0:
        print('Please input valid hmmer score valve.')
        valid = False
    if args.hmmer_e < 0:
        print('Please input valid hmmer e-value limit.')
        valid = False
    return valid


# Annotation group
annotation_parser, annotation_group = register_group('Annotation arguments', [
    {
        'name': 'disable-annotation',
        'default': False,
        'help': 'switching this on will disable annotation.'
    },
    {
        'name': 'species-name',
        'default': 'Test sp.',
        'help': 'species name to use in genbank file.'
    },
    {
        'name': 'use-hmmer',
        'default': False,
        'help': 'if use nhmmer for additional searching for missing PCGs during tblastn, may output low accuracy results, use iif the database is not enough, the results are not ideal and no good way to add profiles.'
    },
    {
        'name': 'hmmer-score',
        'default': 5.0,
        'help': 'if nhmmer is used, results with score higher than this will be selected.'
    },
    {
        'name': 'hmmer-e',
        'default': 0.005,
        'help': 'if nhmmer is used, results with e-value lower than this will be selected.'
    }
], func=annotation_regulator)
