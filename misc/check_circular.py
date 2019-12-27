"""
check_circular.py
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
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import concat_command
except Exception as identifier:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")

from Bio import Seq, SeqIO, SeqRecord
import subprocess


def check_circular(mininum_overlap=8, start_length=300, end_length=300, final_fasta=None):
    try:
        records = list(SeqIO.parse(final_fasta, "fasta"))
    except Exception as i:
        print("Error when loading Fasta file.")
        return None

    dp = os.path.abspath(os.path.join(
        os.path.dirname(__file__), 'dp_circle_check'))

    for record in records:
        seq = str(record.seq)
        f = seq[:start_length]
        r = seq[-end_length:]
        p = subprocess.Popen(dp,
                             stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        result = p.communicate((f + ' ' + r).encode('utf-8')
                               )[0].decode('utf-8').split('\n')[:-1]
        result[0] = [int(x) for x in result[0].split(' ')]
        # TODO: Finish it.
        print(result)
