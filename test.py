from assemble.assemble_wrapper import MEGAHIT
from utility import logger

logger.init('')
megahit = MEGAHIT(
    basedir='../test',
    fq1='../test/test.1.fq',
    fq2='../test/test.2.fq',
    prefix='megatest',
    threads=8
)

kmer_list = [31, 39, 59, 79, 99, 119, 141]
kmer_list = [0, *kmer_list]
megahit.kmax = kmer_list[-1]

megahit.initialize()
megahit.build_lib()

for i in range(1, len(kmer_list)):
    megahit.graph(kmer_list[i - 1], kmer_list[i])
    megahit.assemble(kmer_list[i])
    megahit.filter(kmer_list[i], min_depth=3, min_length=100, max_length=20000)
    if i == len(kmer_list) - 1:
        break
    megahit.local(kmer_list[i], kmer_list[i + 1])
    megahit.iterate(kmer_list[i], kmer_list[i + 1])
