"""
circos.py
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

'''
    An extremely simple tool for ones use python to convert
     a object to a runnable Circos configuration file.

    It's really hacky, but also convinient, the main purpose
    to write this is to make a bridge between python and
    circos, so it may not be very beautiful.

    Equvilant python <=> configuration:

    Python code:
    circos = Circos()

    circos.image.dir = "visualization"
    circos.image.file = "circos.png"
    img = circos.image # Set a reference here to simplify the code.
    img.png = "yes"
    img.svg = "yes"

    circos.show_ticks = "yes"

    Configuration outputted:

    <image>
        dir = visualization
        file = circos.png
        png = yes
        svg = yes
    </image>
    show_ticks = yes
'''


class Circos():

    def __getattribute__(self, name):
        if name != '__dict__':
            if name not in self.__dict__:
                self.__dict__[name] = Circos()
        return object.__getattribute__(self, name)

    def __getitem__(self, index):
        if isinstance(index, tuple):
            key = str(index[0])
            offset = int(index[1])
            if key + '_' * offset in self.__dict__:
                return self.__dict__[key + '_' * offset]
            else:
                return None
        else:
            key = str(index)
            offset = 0
            while key + '_' * offset in self.__dict__:
                offset += 1
            self.__dict__[key + '_' * offset] = Circos()
            return self.__dict__[key + '_' * offset]

    def __enter__(self):
        return self

    def __exit__(self, *_):
        pass


def collapse(obj: Circos):
    dicts = vars(obj)
    for key, value in dicts.items():
        if key.startswith('__'):
            continue
        if isinstance(value, Circos):
            dicts[key] = collapse(value)
    return dicts


def dict2circos(obj: dict, initial: list = None, offset=0):
    if not initial:
        initial = []

    for key, value in obj.items():
        if key.startswith('_'):
            key = key[1:]
        while key.endswith('_'):
            key = key[:-1]

        if isinstance(value, dict):
            initial.append(f'{"    " * offset}<{key}>')
            dict2circos(value, initial, offset+1)
            initial.append(f'{"    " * offset}</{key}>')
        else:
            initial.append(f'{"    " * offset}{key} = {value}')

    return '\n'.join(initial)
