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
    from utility.parser import parse_func, parse_then_call, freeze_main, arg_prop
except ImportError as err:
    sys.exit(f"Unable to import helper module {err}, is the installation of MitoFlex valid?")

from Bio import SeqIO
import subprocess


def check_circular(mininum_length=10000, start_length=500, end_length=500, final_fasta=None):
    try:
        records = list(SeqIO.parse(final_fasta, "fasta"))
    except Exception:
        print("Error when loading Fasta file.")
        return None

    dp = os.path.abspath(os.path.join(
        os.path.dirname(__file__), 'dp_circle_check'))

    finals = []
    for record in records:
        seq = str(record.seq)
        if len(seq) < mininum_length:
            finals.append([-1, record])
            continue
        f = seq[:start_length]
        r = seq[-end_length:]
        p = subprocess.Popen(dp,
                             stdout=subprocess.PIPE, stdin=subprocess.PIPE)
        result = p.communicate((f + ' ' + r).encode('utf-8')
                               )[0].decode('utf-8').split('\n')[:-1]
        result[0] = [int(x) for x in result[0].split(' ')]
        result.append(record)
        finals.append(result)

    return finals


@parse_func
@arg_prop(dest='overlay', default=8, help='the mininum overlay of two sequences')
@arg_prop(dest='length', default=12000, help='the mininum length of sequences, others will be ignored')
@arg_prop(dest='fasta', required=True, default=None, help='input fasta file for sequences being checked')
@arg_prop(dest='start', default=300, help='how many bps will be read at the start of sequence')
@arg_prop(dest='end', default=300, help='how many bps will be read at the end of sequence')
@arg_prop(dest='format', default='std', choices=['std', 'json', 'text'], help='how the script will output the results')
@arg_prop(dest='output', default=None, help='where to write results if there needs one')
def main(args):
    results = check_circular(mininum_length=args.length, start_length=args.start,
                             end_length=args.end, final_fasta=args.fasta)
    if results is None:
        sys.exit(1)

    if args.format == 'std':
        for result in results:
            if result[0] == -1:
                continue
            else:
                if result[0][1] > args.overlay:
                    print(
                        f'Sequence {result[2].name} hit at the position {result[0][0]} with a length of {result[0][1]}.')
                else:
                    print(
                        f'Sequence {result[2].name}\'s overlay is too short and was thought to be linear.')
    else:
        try:
            f = open(args.output, 'w')
        except Exception:
            sys.exit('Errors occured when opening output file.')

        if args.format == 'json':
            # Abusing the generator is fun
            generated = {
                result[2].name: {
                    'start': result[0][0],
                    'length': result[0][1],
                    'overlay': result[1]
                }
                for result in results
                if result[0] != -1 and result[0][1] > args.overlay
            }

            print(generated)

            import json
            json.dump(generated, fp=f)

        elif args.format == 'text':
            # Abusing the generator is really fun...
            [result[0] != -1 and
             result[0][1] > args.overlay and
             print(
                 f'{result[2].name} START: {result[0][0]} LENGTH: {result[0][1]}\nOVERLAY:\n{result[1]}', file=f)
             for result in results]

        f.close()


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
