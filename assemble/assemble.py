"""
assemble.py
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

import os
import sys
try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call, direct_call
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoX valid?")


def assemble(fastq1=None, fastq2=None, result_dir=None, temp_dir=None, work_prefix=None,
             uselist=False, kmin=21, kmax=141, kstep=12, klist=None,
             no_mercy=False, disable_acc=False,
             prune_level=2, prune_depth=2, clean_temp=False,
             threads=8):

    if(uselist):
        kmin = kmax = kstep = None
    else:
        klist = None

    shell_call('megahit', _1=fastq1, _2=fastq2,
               k_min=kmin, k_max=kmax, k_step=kstep, k_list=klist,
               no_mercy=no_mercy, prune_level=prune_level, prune_depth=prune_depth,
               keep_tmp_files=clean_temp, tmp_dir=temp_dir,
               out_dir=result_dir, out_prefix=work_prefix,
               no_hw_accel=disable_acc, num_cpu_threads=threads)

    return os.path.join(result_dir, work_prefix + '.contigs.fa')
