"""
assemble_wrapper.py
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

import os
import sys
from os import path

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, safe_makedirs, timed
    from utility import logger
    from configurations import assemble as a_conf  # Prevent naming confliction
    import psutil
    import uuid
    import subprocess
    from typing import Tuple
except ImportError as err:
    sys.exit(
        f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")


class EmptyGraph(Exception):
    pass


class LibInfo():
    def __init__(self, info):
        self.read_size = int(info[0][0])
        self.read_count = int(info[0][1])
        self.max_len = int(info[2][2])


class ContigInfo():
    def __init__(self, file):
        info = file.readlines()[0].rstrip().split(' ')
        self.count = int(info[0])
        self.bytes = int(info[1])


class MEGAHIT():
    '''
    A more dedicated wrapper to assemble mitogenome sequences with megahit.

    Introduced a filter of depth and length, enabling a much more in-site
    changes of the original pipeline, since MEGAHIT is original designed
    for assembling metagenomes, which are usually of a ultra-low depth in
    sequencing, almost being contradict to the requirements of mitogenome assembly.

    Such a filter can effectively reduce the noise of graph, while retaining
    most of, even improve the quality of assembly, and latter process.
    '''
    basedir = None
    fq1 = None
    fq2 = None
    prefix = None
    threads = None
    use_popcnt = not a_conf.disable_acc
    min_multi = int(a_conf.min_multi)
    no_mercy = bool(a_conf.no_mercy)
    no_local = False
    one_pass = bool(a_conf.one_pass)
    keep_temp = True
    prune_level = 2
    prune_depth = 2
    min_depth = 3
    kmin = 31
    kmax = 141
    additional_kmers = []
    min_length = 200
    max_length = 30000

    def _graph_prefix(self, kmer):
        return path.join(safe_makedirs(path.join(self.temp_dir, f'k{kmer}'), True), str(kmer))

    def _contig_prefix(self, kmer):
        return path.join(self.contig_dir, f'k{kmer}')

    @property
    def MEGAHIT_CORE(self) -> str:
        MEGAHIT_CORES = ["megahit_core_no_hwaccel",
                         "megahit_core_popcnt", "megahit_core"]
        return MEGAHIT_CORES[self.use_popcnt + self.hwaccel]

    @property
    def FAST_FILTER(self) -> str:
        FAST_BIN = path.dirname(__file__)
        return path.join(FAST_BIN, 'fastfilter')

    def __init__(self, **kwargs):
        super().__init__()
        for k, v in kwargs.items():
            self.__dict__[k] = v

    @timed(enabled=False)
    def initialize(self):
        self.basedir = path.abspath(self.basedir)
        self.fq1 = path.abspath(self.fq1)
        if self.fq2:
            self.fq2 = path.abspath(self.fq2)

        self.hwaccel = shell_call("megahit_core checkcpu").rstrip() == '1'
        # Check if POPCNT command is supported
        if self.use_popcnt:
            if shell_call('megahit_core checkpopcnt').rstrip() != '1':
                self.use_popcnt = False
                logger.log(3, "POPCNT is disabled since no features detected.")
            else:
                logger.log(2, f"Using megahit with {'POPCNT' if not self.hwaccel else 'hardware acceleration'} support.")
        else:
            logger.log(2, "POPCNT disabled by argument.")

        if self.one_pass:
            logger.log(3, "Using 1-pass mode.")

        self.result_dir = safe_makedirs(
            path.join(self.basedir, f'{self.prefix}.result'), False)

        if not path.isdir(str(a_conf.external_temp)):
            self.temp_dir = safe_makedirs(
                path.join(self.basedir, f'{self.prefix}.temp'), False)
        else:
            self.temp_dir = safe_makedirs(
                path.join(
                    a_conf.external_temp,
                    str(uuid.uuid4()),
                    f'{self.prefix}.temp'
                ),
                False
            )

        self.read_lib = path.join(self.temp_dir, 'reads.lib')
        self.contig_dir = safe_makedirs(
            path.join(self.temp_dir, 'intermediate_contigs'), False)

        vm = psutil.virtual_memory()
        logger.log(
            1, f"System memory status : {', '.join([f'{k}={v/(1024**2):.2f}MB' for k,v in vm._asdict().items() if type(v) is int])}")
        self.available_memory = int(vm.available * a_conf.max_mem_percent)
        logger.log(2, f'Scheduled {self.available_memory/(1024**2):.2f}MB to use.')

    @timed(enabled=False)
    def build_lib(self):

        # Write reads info
        with open(self.read_lib, 'w') as l:
            fifos = []

            if self.fq1 and self.fq2:
                print(self.fq1, self.fq2, sep=',', file=l)
                fq1, fq2 = (
                    self.fq1 if not self.fq1.endswith('gz') else path.join(self.temp_dir, 'pipe.pe1'),
                    self.fq2 if not self.fq2.endswith('gz') else path.join(self.temp_dir, 'pipe.pe2')
                )

                if self.fq1.endswith('gz'):
                    fifo1 = path.join(self.temp_dir, 'pipe.pe1')
                    os.mkfifo(fifo1)
                    fifos.append(subprocess.Popen(f'gzip -dc {self.fq1} > {fifo1}', shell=True, preexec_fn=os.setsid))

                if self.fq2.endswith('gz'):
                    fifo2 = path.join(self.temp_dir, 'pipe.pe2')
                    os.mkfifo(fifo2)
                    fifos.append(subprocess.Popen(f'gzip -dc {self.fq2} > {fifo2}', shell=True, preexec_fn=os.setsid))

                print('pe', fq1, fq2, file=l)
            else:
                print(self.fq1, file=l)
                fq1 = self.fq1 if not self.fq1.endswith('gz') else path.join(self.temp_dir, 'pipe.se')
                print('se', fq1, file=l)

        logger.log(1, "Converting reads to binary library.")
        shell_call(self.MEGAHIT_CORE, 'buildlib', self.read_lib, self.read_lib)

        if False in (x.wait() == 0 for x in fifos):
            raise RuntimeError("Error occured in reading input fifos")

        with open(self.read_lib + '.lib_info') as ri:
            info = [x.split(' ') for x in ri.readlines()]
            return LibInfo(info)

    @timed(enabled=False)
    def graph(self, current_kmer, next_kmer):
        options = {
            'k': next_kmer,
            'host_mem': self.available_memory,
            'mem_flag': 1,
            'output_prefix': self._graph_prefix(next_kmer),
            'num_cpu_threads': self.threads,
            'need_mercy': not self.no_mercy and current_kmer == self.kmin,
            'kmer_from': current_kmer,
            'useconv': False
        }

        if current_kmer == 0:  # Indicating it's the first graph
            if not self.one_pass:
                logger.log(2, f"Extracting solid (k+1)-mers for k={next_kmer}")
                count_opts = options.copy()
                count_opts['m'] = self.min_multi
                count_opts['read_lib_file'] = self.read_lib
                count_opts.pop('need_mercy')
                count_opts.pop('kmer_from')
                logger.log(0, f"Extract options : {count_opts}")
                shell_call(self.MEGAHIT_CORE, 'count', **count_opts)

        file_size = 0

        if path.exists(self._graph_prefix(next_kmer) + '.edges.0'):
            options['input_prefix'] = self._graph_prefix(next_kmer)
            file_size += path.getsize(self._graph_prefix(next_kmer) + '.edges.0')

        if path.exists(self._contig_prefix(current_kmer) + '.addi.fa'):
            options['addi_contig'] = \
                self._contig_prefix(current_kmer) + '.addi.fa'
            file_size += path.getsize(
                self._contig_prefix(current_kmer) + '.addi.fa')

        if path.exists(self._contig_prefix(current_kmer) + '.local.fa'):
            options['local_contig'] = \
                self._contig_prefix(current_kmer) + '.local.fa'
            file_size += path.getsize(
                self._contig_prefix(current_kmer) + '.addi.fa')

        if path.exists(self._contig_prefix(current_kmer) + '.contigs.fa'):
            options['contig'] = \
                self._contig_prefix(current_kmer) + '.contigs.fa'
            options['bubble'] = \
                self._contig_prefix(current_kmer) + '.bubble_seq.fa'
            file_size += path.getsize(
                self._contig_prefix(current_kmer) + '.contigs.fa')

        if file_size == 0 and current_kmer != 0:
            raise EmptyGraph

        logger.log(2, f'Building graph for k={next_kmer}')
        logger.log(0, f'Build options : {options}')

        shell_call(self.MEGAHIT_CORE, 'seq2sdbg', **options)

        if file_size != 0 and current_kmer != 0 and not self.keep_temp:
            os.system(f"rm -r {path.join(self.temp_dir, f'k{current_kmer}')}")

    @timed(enabled=True)
    def assemble(self, kmer) -> Tuple[ContigInfo, ContigInfo]:
        min_standalone = max(
            min(self.kmax * 3 - 1, int(self.min_length * 1.5)),
            self.min_length)

        options = {
            's': self._graph_prefix(kmer),
            'o': self._contig_prefix(kmer),
            't': self.threads,
            'min_standalone': min_standalone,
            'prune_level': self.prune_level,
            'merge_len': 20,
            'merge_similar': 0.95,
            'cleaning_rounds': 5,
            'disconnect_ratio': 0.1,
            'low_local_ratio': 0.2,
            'min_depth': self.prune_depth,
            'bubble_level': 2,
            'max_tip_len': max(1, self.min_length * 1.5 + 1 - kmer) if kmer * 3 - 1 > self.min_length * 1.5 else -1,
            'careful_bubble': kmer < self.kmax,
            'is_final_round': kmer == self.kmax,
            'output_standalone': self.no_local,
            'useconv': False
        }

        logger.log(2, f'Assembling contigs from SdBG for k = {kmer}')
        logger.log(0, f'Assemble arguments : {options}')

        shell_call(self.MEGAHIT_CORE, 'assemble', **options)
        with open(self._contig_prefix(kmer) + '.contigs.fa.info', 'r') as c, \
                open(self._contig_prefix(kmer) + '.addi.fa.info', 'r') as a:
            return ContigInfo(c), ContigInfo(a)

    @timed(enabled=True)
    def local(self, current_kmer, next_kmer):
        logger.log(2, f'Local assembly for k = {current_kmer}')
        shell_call(self.MEGAHIT_CORE, 'local',
                   c=self._contig_prefix(current_kmer) + '.contigs.fa',
                   l=self.read_lib, t=self.threads,
                   o=self._contig_prefix(current_kmer) + '.local.fa',
                   kmax=next_kmer)

    @timed(enabled=False)
    def iterate(self, current_kmer, next_kmer):
        logger.log(
            2, f'Extracting iterative edges from k = {current_kmer} to {next_kmer}')
        shell_call(self.MEGAHIT_CORE, 'iterate',
                   c=self._contig_prefix(current_kmer) + '.contigs.fa',
                   b=self._contig_prefix(current_kmer) + '.bubble_seq.fa',
                   t=self.threads, s=next_kmer - current_kmer, o=self._graph_prefix(next_kmer),
                   r=self.read_lib + '.bin',
                   k=current_kmer)

    @timed(enabled=False)
    def filter(self, kmer=None,
               min_depth=3, min_length=0, max_length=20000,
               force_filter=False, deny_number=a_conf.filter_keep) -> Tuple[int, int, int]:
        logger.log(2, f'Filtering output contig files of k = {kmer}')

        results = [0, 0, 0]
        if not a_conf.no_filter or force_filter:
            for idx, suffix in enumerate(['.contigs.fa', '.addi.fa', '.bubble_seq.fa']):
                if path.exists(self._contig_prefix(kmer) + suffix):
                    results[idx] = int(
                        shell_call(self.FAST_FILTER,
                                   i=self._contig_prefix(kmer) + suffix,
                                   o=self._contig_prefix(
                                       kmer) + '.filtered' + suffix,
                                   l=f"{min_length},{max_length}",
                                   d=min_depth))

                    if results[idx] <= deny_number and idx == 0:
                        results[idx] = int(shell_call(self.FAST_FILTER,
                                                      i=self._contig_prefix(kmer) + suffix,
                                                      o=self._contig_prefix(kmer) + '.filtered' + suffix,
                                                      l=f"{min_length},{max_length}",
                                                      m=deny_number))

                    shell_call('mv', self._contig_prefix(kmer) + '.filtered' + suffix,
                               self._contig_prefix(kmer) + suffix)

        return tuple(results)

    @timed(enabled=False)
    def finalize(self, kmer):
        self.final_contig = path.join(
            self.result_dir,
            f'k{kmer}.contig.fa'
        )

        shell_call('cat',
                   path.join(self.contig_dir, '*.final.contigs.fa'),
                   self._contig_prefix(kmer) + '.contigs.fa',
                   '>', self.final_contig)

        if not self.keep_temp:
            to_remove = self.temp_dir
            if path.isdir(str(a_conf.external_temp)):
                to_remove = path.join(to_remove, "..")
            to_remove = path.abspath(to_remove)

            os.system(f'rm -r {to_remove}')
