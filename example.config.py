"""
example.config.py
=========

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

# # This is a example configuration, which contains all the argument
# # modifiable in the MitoFlex.
# # The configuration file doesn't have a mechanism to check out which
# # one is the really 'valid' data, configuration will always override
# # command-line values, so be catious in checking values. The only one
# # way NOT to overwrite arguments passed is to comment out the values
# # you don't want to.
# # Also, this configuration used a trick to make python import and
# # execute it as an module, thus you may add your config with some
# # code to make it even more flexible...

# # Global configuration
# command = 'all'
#
# # Universal configuration
# workname = ''
# threads = 8
#
# # Fastq configuration
# fastq1 = ''
# fastq2 = ''
# fastq_alter_format = False
# fastq_read_length = 150
# fq_size = 5
#
# # Fasta configuration
# fastafile = ''
#
# # Filter configuration
# adapter1 = ''
# adapter2 = ''
# cleanq1 = ''
# cleanq2 = ''
# deduplication = False
# adapter_mismatch = 3
# adapter_length = 15
# Ns_valve = 10
# quality_valve = 55
# percentage_valve = 0.2
# keep_region = ''
#
# #Assembly configuration
# insert_size = 150
# kmer_min = 31
# kmer_max = 63
#
# # Search mitochondrial gene configuration
# filter_taxa = False
# min_abundance = 10
# required_taxa = 'Platyhelminthes'
# taxa_tolerance = 0
#
# # Search and annotation mitochondrial gene configuration
# genetic_code = 9
# clade = 'Platyhelminthes-flatworms'
#
# #Annotation group
# disable_annotation = True
# species_name = 'Test sp.'
