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
from typing import List

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.parser import parse_func, parse_then_call, freeze_main, arg_prop
    from misc import libfastmathcal
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

except ImportError as err:
    sys.exit(f"Unable to import helper module {err}, is the installation of MitoFlex valid?")


def check_circular(mininum_length=10000, start_length=500, end_length=500, overlaps=50, final_seqs: List[SeqRecord] = None):
    if final_seqs is None:
        return

    for record in final_seqs:
        seq = str(record.seq)
        if len(seq) < mininum_length:
            yield (None, record)
            continue
        f = seq[:start_length]
        r = seq[-end_length:]
        f_start, f_end, ali_length = libfastmathcal.seq_overlap(f, r)
        if ali_length < overlaps:
            yield (None, record)
        else:
            yield ((f_start, f_end, ali_length), record)


@parse_func
@arg_prop(dest='overlay', default=8, help='the mininum overlay of two sequences')
@arg_prop(dest='length', default=12000, help='the mininum length of sequences, others will be ignored')
@arg_prop(dest='fasta', required=True, default=None, help='input fasta file for sequences being checked')
@arg_prop(dest='start', default=300, help='how many bps will be read at the start of sequence')
@arg_prop(dest='end', default=300, help='how many bps will be read at the end of sequence')
@arg_prop(dest='output', default=None, help='where to write results if there needs one')
def main(args):
    results = {x[1].id: x[0] for x in check_circular(args.length, args.start, args.end, args.length, SeqIO.parse(args.fasta, 'fasta'))}
    with open(args.output, 'w') as f:
        import json
        json.dump(results, f)


desc = '''
check_circular.py

Description
    This is a script written to check if a input of sequences could
    be a circle or just linear.
    Input fasta file, print found sequences, or you can set output to
    json or something else.
    This uses an algorithm of O(n^2) time complexity, and should be
    fast enough with any inputs with search regions from 300 to 3k bps.
'''

# Program entry starts at here
if __name__ == '__main__':
    parser = freeze_main(prog='check_circular.py', desc=desc)
    parse_then_call(parser)
