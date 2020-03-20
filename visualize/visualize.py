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
    from utility.bio import circos
    from visualize import circos_config
    from utility import logger
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")


def visualize(fasta_file=None, fastq1=None, fastq2=None, pos_json=None,
              prefix=None, basedir=None, threads=8, circular=False):
    logger.log(2, 'Entering visualize module.')
    # Validate the paths
    fasta_file = path.abspath(fasta_file)
    fastq1 = path.abspath(fastq1)
    fastq2 = path.abspath(fastq2)
    basedir = path.abspath(basedir)
    pos_json = path.abspath(pos_json)

    fa_copy = path.join(basedir, f'{prefix}.fasta')
    list_conv = []
    counter = 1

    # Rename to a easier form
    index_list = {}
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        index_list[seq.id] = f'mt{counter}'
        seq.id_old = seq.id
        seq.id = f'mt{counter}'
        seq.description = ''
        list_conv.append(seq)
        counter += 1
    SeqIO.write(list_conv, fa_copy, 'fasta')

    with open(pos_json, 'r') as f:
        poses = json.load(f)

    # Gene name files
    logger.log(1, 'Generating gene name and feature files.')
    gene_name_file = path.join(basedir, f'{prefix}.gene.txt')
    with open(gene_name_file, 'w') as gn_f:
        for key, value in poses.items():
            start, end, gene_type, strand = value
            strand_conv = index_list[strand]
            print(strand_conv, start, end, key, sep='\t', file=gn_f)

    # Gene feature files
    gene_feature_file = path.join(basedir, f'{prefix}.features.txt')
    with open(gene_feature_file, 'w') as gf_f:
        for key, value in poses.items():
            start, end, gene_type, strand = value
            strand_conv = index_list[strand]
            print(strand_conv, start, start,
                  'fill_color=black,r0=0.965r,r1=1r', file=gf_f, sep='\t')
            print(strand_conv, start, end,
                  f'fill_color={circos_config.fill_colors[int(gene_type)]},r0=0.965r,r1=1r',
                  file=gf_f, sep='\t')
            print(strand_conv, end, end,
                  'fill_color=black,r0=0.965r,r1=1r', file=gf_f, sep='\t')

    logger.log(1, 'Generating depth files.')
    # Using check_output directly because being too lazy to remove decoder
    from subprocess import check_output

    shell_call('bwa index', fa_copy)
    bam_file = path.join(basedir, f'{prefix}.bam')
    check_output(
        f'bwa mem -t {threads} {fa_copy} {fastq1} {fastq2} |samtools view -bS -q 30 -h -o {bam_file} -', shell=True)
    bam_sorted_file = path.join(basedir, f'{prefix}.sorted.bam')
    check_output(f'samtools sort -o {bam_sorted_file} {bam_file}', shell=True)
    gene_depth_file = path.join(basedir, f'{prefix}.dep')
    check_output(
        f'samtools depth -aa {bam_sorted_file} > {gene_depth_file}', shell=True)

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
    # Reusing conv-list here, as it's not deleted in the scope
    gc_content_file = path.join(basedir, f'{prefix}.gc.txt')
    with open(gc_content_file, 'w') as gc_f:
        for seq in list_conv:
            # Stepping 50 to walk through
            for s in range(0, len(seq), 50):
                seq_slice = seq[s:s+50]
                gc_num = sum(x == 'G' or x == 'C' for x in seq_slice)
                gc_per = gc_num / len(seq_slice)
                print(seq.id, s, s+len(seq_slice), gc_per, file=gc_f)

    # Karyotype
    logger.log(1, 'Generating chr files.')
    karyotype_file = path.join(basedir, f'{prefix}.karyotype.txt')
    with open(karyotype_file, 'w') as ky_f:
        for seq in list_conv:
            chr_name = seq.id.replace('mt', 'chr')
            print(
                f'{chr_name} - {seq.id}\t{seq.id_old}\t0\t{len(seq)}\tgrey', file=ky_f)

    # Plus generation
    logger.log(1, 'Generating plus.')
    plus_file = path.join(basedir, f'{prefix}.plus.txt')
    with open(plus_file, 'w') as p_f:
        print('mt1\t0\t300\t+\tr0=1r-150p,r1=1r-100p', file=p_f)

    # Giving the values
    logger.log(1, 'Generating circos config file.')
    generated_config = circos_config.circos_conf
    generated_config.ideogram.spacing._break = "0.5r" if not circular else "0.01r"
    generated_config.image.dir = basedir
    generated_config.karyotype = karyotype_file
    generated_config.plots['plot', 0].file = gene_name_file
    generated_config.plots['plot', 1].file = plus_file
    generated_config.plots['plot', 2].file = gc_content_file
    with generated_config.plots['plot', 3] as depth_plot:
        depth_plot.file = circos_depth_file
        depth_plot.max = max_gene_depth
        depth_plot.rules['rule', 0
                         ].condition = f'var(value) > {int(max_gene_depth*0.9)}'
        depth_plot.rules['rule', 1
                         ].condition = f'var(value) < {int(max_gene_depth*0.1)}'

    generated_config.highlights['highlight', 0].file = gene_feature_file

    # Writing to final
    # I guess it would be better to use a f-string formatted cfg, but
    # well this is fine.
    cfg_dict = circos.collapse(generated_config)
    cfg_file = path.join(basedir, 'circos.conf')
    with open(cfg_file, 'w') as cfg_f:
        cfg_f.write('<<include etc/colors_fonts_patterns.conf>>\n')
        cfg_f.write(circos.dict2circos(cfg_dict) + '\n')
        cfg_f.write('<<include etc/housekeeping.conf>>')

    logger.log(1, 'Running Circos.')
    check_output('circos', shell=True, cwd=basedir)
    return path.join(basedir, 'Circos.png'), path.join(basedir, 'Circos.svg')
