"""
configurations.py
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


'''
This file is for storing things that can be configurated, but too
static and less used, like some experimental methods, strange modules
and functions that behaves badly in most of the time, and things that just
too fit to be changed.
But this may work well in certain situations, so you may will want to modify
the value inside this configuration if needed.
'''

# It's not acutally circos configuration file, but used to make the file structure clean
from utility.bio.circos import Circos
annotation = Circos()
findmitoscaf = Circos()
assemble = Circos()
filter_rawdata = Circos()
visualize = Circos()

# Filter

# Enable clean data compression?
# Enabling the data compression will make the filter to automatically
# compress the output data, but will also have a huge performance impact
# on the process time (5min -> 30min, even 1h)
filter_rawdata.compress_output_in_all = False

# Assemble

# Setting the max memory percent for MEGAHIT to use .
# MEGAHIT will NOT allocate all the memory at once, but will only allocate
# this percentage at most.

assemble.max_mem_percent = 0.9

# Setting the min multiplicity of k-mer input.
# Low value works well if the input reads are of low depth and coverage,
# but obviously mitogenome shouldn't.
assemble.min_multi = 3

# Do not allow mercy edges?
# Disable this will allow mercy edges in graph. Works well if input reads
# are of low depth and coverage, but will increase the noise and complexity
# of the graph, making this option rather useless when assembling mitogenome/
assemble.no_mercy = True

# K-mer 1-pass mode?
# Enabling 1-pass mode will make the megahit performs better in situations of
# ultra-low gene coverage, also makes it useless in MitoFlex.
assemble.one_pass = False

# Disable the POPCNT support of megahit?
# May have a performance impact. But it may makes the assembler work on old
# CPUs, only use as a last resort.
assemble.disable_acc = False

# Disable Fast-Filter for controlling MEGAHIT's assemble results?
# Use for concentrating graph to sequences more likely to be mitogenome's sequences.
# Disabling it will make MEGAHIT work as original, which may have bad results on
# mitogenome sequences.
assemble.no_filter = False

# The minimum contigs to keep when assembling contigs.
# Set 0 to always filter, if output contigs is less than this number, filter will no
# be applied.
assemble.filter_keep = 5000

# What contig matching the conditions will be retained in assembly?
# Use to specify how you want to control your contig output, more strict
# values will lead to more accurate but less results, and vise versa.

# The minimum length for contig to output in the iteration.
assemble.min_length = 200

# The maximum length for contig to output in the iteration.
assemble.max_length = 20000

# Show warning of multi value outputted by SOAP-Wrapper?
assemble.show_from_soap = True

# The max number of threads when scaffolding.
# This argument is for limiting the thread used in the scaffolding, since
# it's not completely thread safe when mapping sequences, specifying a higher
# number may increase the speed of mapping and scaffolding, while have a higher
# risk to halt the progress, or output few even no scaffolds.
assemble.max_thread_scaf = 16

# If use external (and fixed) temp directory and where.
# None means do not use, and if a valid directory is specified, assembler will use
# that instead. This is for somewhere like tmpfs can then be utilized, to speed
# up the IO.
# The external buffer needs about 20G or more space to store the temporal data.
# And IO rate may have little improvement if :
# 1. You are assembling a very large genome.
# 2. Comparing to your disk speed, your memcpy speed is relatively slow.
#    This can occur if you have a RAID or SSD but a outdated RAM, like DDR3.
assemble.external_temp = None

# Findmitoscaf

# Default clade of data.
# Change this to the clades you frequently studies to avoid explicitly
# specify the clade every time.
findmitoscaf.default_clade = "Platyhelminthes"

# What partial do MitoFlex treat the gene of a single sequence as full gene?
# In most circumstances, genes are just full, but in case that your clade is rare,
# and have no precise HMM models, lower this may increase the accuracy of the contig
# result.

findmitoscaf.full_ratio = 0.95

# What partial do MitoFlex treat the alignment of a single sequence as a valid gene?
# MitoFlex will think of some align in center as full if align is longer than this
# ratio, not considering if the thing meets the full_ratio. Adding this may increase
# the sensitivity, while lowering the accuracy.

findmitoscaf.min_valid_ratio = 0.3

# Should another findmitoscaf run to be launched after the merging.
# Since some sequences will be conflicted, or resulted to be have no gene at all after
# the merge method, an additional check can improve the result quality, though at some
# cost of losing already found genes. If you are feeling like the risk of potentially
# overlapping is not that bad, and losing gene is quite unbearable, please turn this off.

findmitoscaf.additional_check = True

# Should findmitoscaf sequence be merged broke into two parts?
# It's specialized for annotating genes that are happened to be splited at the start,
# or the end of the linear mitogenome.
# Adding this will make MitoFlex to add an additional sequence that merges the start
# and the end of the sequence, if the picked mitoscaf is only 1 sequence and circular.
# Extra sequence will be marked to be disposed at the end of annotation, so there will
# be no chance in to result.
findmitoscaf.split_two = False

# Annoation

# Enable gene search relocation?
# Enabling relocation will make Mitoflex tries to find the accurate start and end
# codon of each PCG, may not working well on certain situations.
annotation.reloc_genes = False

# Enable genome trimming if it's a circular sequence?
# Disabling the genome trimming will also disable the program's ability of detecting
# circular genome sequences, where the output will always be linear.
annotation.trim_circular = True

# Enable genome redirecting if many of the PCGs are found in another strand
# This will check if the genome is have more than 50% of the gene in positive
# strand, and will reverse it before annotation if not.
annotation.redirection = False

# The ratio of two 'valid' overlapping gene
# MitoFlex will remove overlapped tblastn results if the overlapped region is
# longer than the ratio*min(seq1, seq2). Increasing this will have more results
# with lower accuracy, and vise versa.
annotation.overlap_ratio = 0.2

# The fill color of genes in visualize method
visualize.color_cds = '141,211,199'
visualize.color_trna = '251,128,114'
visualize.color_rrna = '253,192,134'
