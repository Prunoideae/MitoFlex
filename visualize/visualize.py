"""
visualize.py
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

import os
import json
import sys
import shutil
from os import path

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call
    from Bio import SeqIO
    from visualize import circos_config
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")


def visualize(fasta_file=None, fastq1=None, fastq2=None, pos_json=None,
              prefix=None, basedir=None, threads=8):
    # Validate the paths
    fasta_file = path.abspath(fasta_file)
    fastq1 = path.abspath(fastq1)
    fastq2 = path.abspath(fastq2)
    basedir = path.abspath(basedir)
    pos_json = path.abspath(pos_json)

    fa_copy = path.join(basedir, f'{prefix}.{path.splitext(fasta_file)[1]}')
    list_conv = []
    counter = 1

    # TODO test after the dataset is done
    # Rename to a easier form
    index_list = {}
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        index_list[seq.id] = f'mt{counter}'
        seq.id = f'mt{counter}'
        seq.description = ''
        list_conv.append(seq)
        counter += 1
    SeqIO.write(list_conv, fa_copy, 'fasta')

    with open(pos_json, 'r') as f:
        poses = json.load(pos_json)

    redirected_pos = {index_list[key]: value for key, value in poses.items()}

    shell_call('bwa index', fa_copy)
    bam_file = path.join(basedir, f'{prefix}.bam')
    shell_call(f'bwa mem -t {threads} "{fasta_file}", "{fastq1}", "{fastq2}" |',
               'samtools view -bS', q=30, h=True, o=bam_file)
    bam_sorted_file = path.join(basedir, f'{prefix}.sorted.bam')
    shell_call('samtools sort', bam_file, o=bam_sorted_file)
    gene_depth_file = path.join(f'{prefix}.dep')
    shell_call('samtools depth -aa', bam_sorted_file, '>', gene_depth_file)

    # Calculate the things
    circos_depth_file = path.join(basedir, f'{prefix}.depth.txt')
    max_gene_depth = 0
    with open(gene_depth_file, 'r') as gdf, open(circos_depth_file, 'w') as cdf:
        for line in gdf:
            content = str(line).rstrip().split()
            print(' '.join([content[0], content[1],
                            content[1], content[2]]), file=cdf)
            if int(content[2]) > max_gene_depth:
                max_gene_depth = int(content[2])

    # GC content
    # Reuse conv-list here, as it's not deleted in the scope
    for seq in list_conv:
        # Stepping 50 to walk through
        pass

    # Giving the values
    generated_config = circos_config.circos_conf
    generated_config.image.dir = 'Output directory here'
    generated_config.karyotype = 'Karyotype file here'
    generated_config.plots['plot', 0].file = 'Gene name and position file here'
    generated_config.plots['plot', 1].file = 'Plus file here'
    generated_config.plots['plot', 2].file = 'GC content file here'
    with generated_config.plots['plot', 3] as depth_plot:
        depth_plot.file = 'Depth file here'
        depth_plot.rules['rule', 0].condition = 'var(value) > {}'
        depth_plot.rules['rule', 1].condition = 'var(value) < {}'

    generated_config.highlights['highlight', 0].file = "Feature file here"
