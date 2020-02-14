"""
infernal.py
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
import sys

try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    #import wuss
    from utility.bio import wuss
except Exception as iden:
    sys.exit("Unable to import helper module, is the installation of MitoFlex valid?")


'''
    A simple parser to parse output from Infernal.
    This parser only process one-line results, as the intention to write
    it is dedicated.
'''


class Infernal():

    class Result():
        def __init__(self, data: str):
            lines = data.split('\n')[0:-1]
            ids = lines[0].split()
            self.sequence = ids[0]
            paras = lines[3].split()
            self.rank = paras[0].translate({ord(i): None for i in '()'})
            self.e_value = float(paras[2])
            self.score = float(paras[3])
            self.bias = float(paras[4])
            self.mdlfrom = int(paras[6])
            self.mdlto = int(paras[7])
            self.seqfrom = int(paras[9])
            self.seqto = int(paras[10])
            self.plus = paras[11] == '+'
            self.acc = float(paras[13])
            self.full = paras[14] == 'no'
            self.gc = float(paras[15])
            self.alignment = None

            seq = lines[8].split(maxsplit=2)[2].rsplit(maxsplit=1)[0]
            fold = lines[5].split()[0]

            fold, seq = wuss.align_fold(fold, seq)

            sing = wuss.seq2single(seq)
            self.qual = lines[9].split()[0]
            self.alignment = wuss.GenericLoop(fold, sing)

        def __repr__(self):
            return '\n'.join([f'{key}={value}' for key, value in vars(self).items()])

        def __str__(self):
            return '\n'.join([f'{key}={value}' for key, value in vars(self).items()])

    def __init__(self, file):

        alignments = []

        # Read file into result
        try:
            with open(file, 'r') as inf:
                stage = 0
                for line in inf:
                    if line != '\n':
                        if line.startswith('#'):
                            continue
                        elif line.startswith('Hit alignments'):
                            stage = 1
                            continue
                        elif '[No hits detected that satisfy reporting thresholds]' in line:
                            continue
                        elif line.startswith('Internal CM pipeline statistics summary'):
                            stage = 2
                            continue

                        if stage == 1:
                            alignments.append(line)
        except Exception:
            raise IOError("Cannot read infernal file!")

        alignments = ''.join(alignments).split('>> ')[1:]
        self.alignments = [Infernal.Result(
            data) for data in alignments if '[No hits detected that satisfy reporting thresholds]' not in data]


class Queries():

    class Query():
        def __init__(self, data: str):
            splited = data.split()
            self.rank = int(splited[0].translate({ord(i): None for i in '()'}))
            self.e_value = float(splited[2])
            self.score = float(splited[3])
            self.bias = float(splited[4])
            self.sequence = str(splited[5])
            self.seqfrom = int(splited[6])
            self.seqto = int(splited[7])
            self.plus = splited[8] == '+'
            self.full = splited[10] == 'no'
            self.gc = float(splited[11])

    def __init__(self, data: str):
        queries = []
        try:
            with open(data, 'r') as f:
                stage = 0
                for line in f:
                    if line != '\n':
                        if line.startswith('#'):
                            continue
                        elif line.startswith('Query'):
                            stage = 1
                            continue
                        elif line.startswith('Hit alignments'):
                            break

                        if stage == 1:
                            queries.append(line)

        except Exception:
            raise IOError("Cannot read infernal file!")

        queries = queries[3:]
        queries = ''.join(queries).split('\n')[:-1]

        self.queries = [Queries.Query(
            x) for x in queries if '[No hits detected that satisfy reporting thresholds]' not in x]
