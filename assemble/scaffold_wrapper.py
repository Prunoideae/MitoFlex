import os
import sys
from os import path

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call
    from misc.check_circular import check_circular
    from utility import logger
    from Bio import SeqIO
    from configurations import assemble as a_conf
except ImportError as err:
    sys.exit(
        f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

bin_dir = path.dirname(__file__)
soap_fusion = path.join(bin_dir, 'SOAPdenovo-fusion')
soap_127 = path.join(bin_dir, 'SOAPdenovo-127mer')


class SOAP():

    def __init__(self, fq1, fq2, contigs, read_length, insert_size, basedir, threads, prefix, final_kmer):
        self.fq1 = path.abspath(fq1)
        self.fq2 = path.abspath(fq2)
        self.contigs = path.abspath(contigs)
        self.insert_size = insert_size
        self.lib_file = None
        self.basedir = path.join(path.abspath(basedir), f"{prefix}.scaf")
        self.read_length = read_length
        self.threads = min(threads, a_conf.max_thread_scaf)
        self.final_kmer = final_kmer
        os.mkdir(self.basedir)

    def lib(self):
        avg_ins = self.insert_size
        self.lib_file = path.join(self.basedir, "soaplib.txt")
        with open(self.lib_file, 'w') as f:
            f.write(
                f'max_rd_len={self.read_length}\n'
                '[LIB]\n'
                f'avg_ins={avg_ins}\n'
                'reverse_seq=0\n'
                'asm_flags=3\n'
                'rank=1\n'
                'pair_num_cutoff=3\n'
                'map_len=32\n'
                f'q1={self.fq1}\n' + f"q2={self.fq2}" if self.fq2 != None else "")

    def scaf(self) -> str:
        if self.lib_file == None:
            raise RuntimeError("Lib was not build before scaffolding!")

        kmer = int(self.read_length / 2)
        prefix = path.join(self.basedir, f'k{kmer}')

        # Prepare
        logger.log(2, "Constructing graph for SOAPdenovo-127.")
        shell_call(soap_fusion, D=True, s=self.lib_file,
                   p=self.threads, K=kmer, g=prefix, c=self.contigs)

        # Map
        logger.log(2, "Mapping sequences.")
        shell_call(soap_127, 'map', s=self.lib_file,
                   p=self.threads, g=prefix)

        # Scaff
        logger.log(2, "Scaffolding.")
        shell_call(soap_127, 'scaff', p=self.threads, g=prefix)

        # Convert
        logger.log(2, "Converting output scaffolds back.")
        scaf2mega(prefix + '.scafSeq',
                  path.join(path.dirname(self.contigs), 'scaf.fa'),
                  overlay=kmer)
        return path.join(path.dirname(self.contigs), 'scaf.fa')


def scaf2mega(i, o, overlay):
    if a_conf.show_from_soap:
        logger.log(
            3, "NOTICE: due to the limit of SOAPdenovo-fusion and 127mer, scaffolds' depths are not correctly calculated.")
        logger.log(
            3, "To avoid the later process to unwisely filter out scaffolds, these sequences are always tolerated!")
        logger.log(
            3, "But don't worry, if you have a correct depth filter setup, output scaffolds should always be safe enough.")
        logger.log(
            3, "You can disable this message in the configurations.py if you have already knew this.")

    with open(o, 'w') as f:
        for info, seq in check_circular(1000, overlay * 4, overlay * 4, overlaps=overlay, final_seqs=SeqIO.parse(i, 'fasta')):
            flag = 1 if info is None else 3
            seq.description = f'flag={flag} multi=32767 len={len(seq)}'
            SeqIO.write(seq, f, 'fasta')
