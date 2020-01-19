"""
wuss.py
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


from itertools import groupby


'''
A simple parser for parsing Washington University
Secondary Structure (WUSS) parsing.
How did you have this name as the format...?
I'm confused.
'''


# Conversion

def seq2single(sequence: str):
    return [Single(x, []) for x in sequence]


# Foundamental structures

class Single():
    def __init__(self, base, parent: list):
        self.base = base
        self.parent = parent

    def __repr__(self):
        return self.base


class Sequence():
    def __init__(self, sequence: list):
        self.sequence = sequence

    def push(self, base: Single):
        self.sequence.append(base)

    def __repr__(self):
        return ''.join([single.base for single in self.sequence])


class Sets():
    def __init__(self, bases: set):
        self.bases = bases

    def insert(self, base: Single):
        self.bases.add(base)

    def __repr__(self):
        return f'({",".join([single.base for single in self.bases])})'


class Paired():
    def __init__(self, left: list, right: list):
        self.left = left
        self.right = right

    def insert(self, l: Single, r: Single):
        self.left.insert(0, l)
        self.right.append(r)

    def __repr__(self):
        return f'L:{"".join([single.base for single in self.left])} R:{"".join([single.base for single in self.right])}'


# Components

class Hairpin(Sequence):
    def __init__(self, sequence):
        super().__init__(sequence)


class Stem(Paired):
    def __init__(self, left, right):
        super().__init__(left, right)


class InteriorLoop(Sets):
    def __init__(self, bases: set):
        super().__init__(bases)


class MultiBranchLoop(Sets):
    def __init__(self, bases: set):
        super().__init__(bases)


# Partions of Secondary Structures

class HairpinLoop():

    '''
    Any loops that enclosed by <>, which indicates it a hairpin loop.
    '''

    def __init__(self, foldseq: str, sequence: list):
        if len(foldseq) != len(sequence):
            raise RuntimeError(
                "Fold sequence should be as long as the base sequence!")
        self.sequence = sequence
        self.loop = InteriorLoop(set())
        self.hairpin = Hairpin([])
        self.stem = Stem([], [])
        stack = []
        for idx, cha in enumerate(foldseq):
            sequence[idx].parent.append(self)
            if cha == '_':
                sequence[idx].parent.append(self.hairpin)
                self.hairpin.push(sequence[idx])
            elif cha == '<':
                sequence[idx].parent.append(self.stem)
                stack.append(sequence[idx])
            elif cha == '>':
                sequence[idx].parent.append(self.hairpin)
                ident = stack.pop()
                self.stem.insert(ident, sequence[idx])
            elif cha == '-':
                sequence[idx].parent.append(self.loop)
                self.loop.insert(sequence[idx])
        son_level = sequence[0].parent.index(self)+1
        translated = [base.parent[son_level] for base in sequence]
        self.components = [x[0] for x in groupby(translated)]


class MultiLoop():
    '''
    Any loops that is enclosed by (), indicating the loop has serveral hairpins or others.
    '''

    def __init__(self, fold: str, sequence: list):
        if len(fold) != len(sequence):
            raise RuntimeError(
                "Fold sequence should be as long as the base sequence!")
        self.sequence = sequence
        self.stem = Stem([], [])
        self.multi = MultiBranchLoop(set())
        self.interior = InteriorLoop(set())

        stack = []

        for idx, cha in enumerate(fold):
            sequence[idx].parent.append(self)

            if cha == '(':
                stack.append(('(', sequence[idx]))
            elif cha == ')':
                ident = stack.pop()[1]
                self.stem.insert(ident, sequence[idx])
                sequence[idx].parent.append(self.stem)
                ident.parent.append(self.stem)
            elif cha == '<':
                stack.append(('<', idx))
            elif cha == '>':
                pos = stack.pop()[1]
                if '<' not in [x[0] for x in stack]:
                    HairpinLoop(fold[pos:idx+1], sequence[pos:idx+1])
            elif cha == ',':
                self.multi.insert(sequence[idx])
                sequence[idx].parent.append(self.multi)
            elif cha == '-':
                if '<' not in [x[0] for x in stack]:
                    self.interior.insert(sequence[idx])
                    sequence[idx].parent.append(self.interior)

        son_level = sequence[0].parent.index(self)+1
        translated = [base.parent[son_level] for base in sequence]
        self.components = [x[0] for x in groupby(translated)]


class ComplexLoop():
    '''
    Any loop that is enclosed by [], indicating it a higher level than MultiLoop
    '''

    def __init__(self, fold: str, sequence: list):
        if len(fold) != len(sequence):
            raise RuntimeError(
                "Fold sequence should be as long as the base sequence!")
        self.sequence = sequence
        self.multi = MultiBranchLoop(set())
        self.stem = Stem([], [])
        self.interior = InteriorLoop(set())
        self.mismatch = Sets(set())

        stack = []

        for idx, cha in enumerate(fold):
            sequence[idx].parent.append(self)

            if cha == '[':
                stack.append(('[', sequence[idx]))
            elif cha == ']':
                ident = stack.pop()[1]
                self.stem.insert(ident, sequence[idx])
                ident.parent.append(self.stem)
                sequence[idx].parent.append(self.stem)
            elif cha == '(':
                stack.append(('(', idx))
            elif cha == ')':
                ident = stack.pop()[1]
                if '(' not in [x[0] for x in stack]:
                    MultiLoop(fold[ident: idx+1], sequence[ident:idx+1])
            elif cha == '<':
                if '(' not in [x[0] for x in stack]:
                    stack.append(('<', idx))
            elif cha == '>':
                if '<' in [x[0] for x in stack]:
                    ident = stack.pop()[1]
                    if '<' not in [x[0] for x in stack]:
                        HairpinLoop(fold[ident:idx+1], sequence[ident:idx+1])
            elif cha == ',':
                if '(' not in [x[0] for x in stack]:
                    self.multi.insert(sequence[idx])
                    sequence[idx].parent.append(self.multi)
            elif cha == '-':
                if not(False in [x[0] == '[' for x in stack]):
                    self.interior.insert(sequence[idx])
                    sequence[idx].parent.append(self.interior)
            elif cha == ':':
                self.mismatch.insert(sequence[idx])
                sequence[idx].parent.append(self.mismatch)

        son_level = sequence[0].parent.index(self)+1
        translated = [base.parent[son_level] for base in sequence]
        self.components = [x[0] for x in groupby(translated)]


# General structure

class GenericLoop():

    '''
    Any loops that is enclosed by {}, means that the loop is on a even higher level.
    Also, it can parse loops that of any lower level, since it's the highest one of
    all, the name 'GenericLoop' is called.
    '''

    def __init__(self, fold: str, sequence: list):
        if len(fold) != len(sequence):
            raise RuntimeError(
                "Fold sequence should be as long as the base sequence!")
        self.sequence = sequence
        self.multi = MultiBranchLoop(set())
        self.stem = Stem([], [])
        self.interior = InteriorLoop(set())
        self.mismatch = Sets(set())

        stack = []

        for idx, cha in enumerate(fold):
            sequence[idx].parent.append(self)

            if cha == '{':
                stack.append(('{', sequence[idx]))
            elif cha == '}':
                ident = stack.pop()[1]
                self.stem.insert(ident, sequence[idx])
                ident.parent.append(self.stem)
                sequence[idx].parent.append(self.stem)
            elif cha == '[':
                stack.append(('[', idx))
            elif cha == ']':
                ident = stack.pop()[1]
                if '[' not in [x[0] for x in stack]:
                    ComplexLoop(fold[ident:idx+1], sequence[ident:idx+1])
            elif cha == '(':
                if '[' not in [x[0] for x in stack]:
                    stack.append(('(', idx))
            elif cha == ')':
                if '(' in [x[0] for x in stack]:
                    ident = stack.pop()[1]
                    if '(' not in [x[0] for x in stack]:
                        MultiLoop(fold[ident: idx+1], sequence[ident:idx+1])
            elif cha == '<':
                if not(False in [x[0] == '{' or x[0] == '<' for x in stack]):
                    stack.append(('<', idx))
            elif cha == '>':
                if '<' in [x[0] for x in stack]:
                    ident = stack.pop()[1]
                    if '<' not in [x[0] for x in stack]:
                        HairpinLoop(fold[ident:idx+1], sequence[ident:idx+1])
            elif cha == ',':
                if not(False in [x[0] == '{' for x in stack]):
                    self.multi.insert(sequence[idx])
                    sequence[idx].parent.append(self.multi)
            elif cha == '-':
                if not(False in [x[0] == '{' for x in stack]):
                    self.interior.insert(sequence[idx])
                    sequence[idx].parent.append(self.interior)
            elif cha == ':':
                self.mismatch.insert(sequence[idx])
                sequence[idx].parent.append(self.mismatch)

        son_level = sequence[0].parent.index(self)+1
        translated = [base.parent[son_level] for base in sequence]
        self.components = [x[0] for x in groupby(translated)]
