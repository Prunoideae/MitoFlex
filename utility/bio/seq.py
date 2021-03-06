"""
seq.py
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

# Helper class to make code struct clean


# compile a dict to a description
def compile_seq(traits: dict, equ='=', sep='\t'):
    return sep.join(equ.join((key, str(value))) for key, value in traits.items())


# decompile a description to a dict
def decompile(input_seq: str, equ='=', sep='\t'):
    return {trait.split(equ)[0]: trait.split(equ)[1]
            for trait in input_seq.split(sep)
            if equ in trait}
