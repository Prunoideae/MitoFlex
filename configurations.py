"""
configurations.py
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


'''
This file is for storing things that can be configurated, but too
static and less used, like some experimental methods, strange modules
and process that behaves badly in most of the time.
'''

# It's not acutally circos configuration file, but used to make the code structure clean
from utility.bio.circos import Circos
annotation = Circos()
findmitoscaf = Circos()
assemble = Circos()
filter_rawdata = Circos()
visualize = Circos()

# Filter

# Enable clean data compression?
filter_rawdata.compress_output_in_all = False

# Annoation

# Enable gene search relocation?
annotation.reloc_genes = False
