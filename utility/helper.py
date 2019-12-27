"""
helper.py
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

import subprocess
import sys
import os

# cmd runner


def shell_call(*args, **kwargs):
    '''
    Concatenacte all the args and kwargs into a shell command string, then run it.\n
    The input arguments are not limited, but rules should be told.
    1. Positional arguments will always comes before the keyword arguments.
    2. This func will always add '--' before the keyword, maybe other prefixes will be
    acceptable but not now.
    shell_call('python', 'fun.py', foo='bar') -> 'python fun.py --foo bar'
    shell_call('python', 'foo.py', bar='lorem ipsum') -> 'python foo.py --bar lorem ipsum'
    shell_call('python', 'bar.py', wow_fun='method') -> 'python bar.py --wow-fun method'
    '''
    command = concat_command(*args, **kwargs)
    direct_call(command)


def concat_command(*args, **kwargs):
    args = [str(x) for x in args]
    kwargs = {x[1:]if x.startswith('_') else x: kwargs[x]
              for x in kwargs if kwargs[x] is not None}
    kwargs = {str(x).replace('_', '-'): kwargs[x] for x in kwargs}
    command = ' '.join(args)
    for arg in kwargs:
        if kwargs[arg] is None:
            continue
        if type(kwargs[arg]) is list:
            command += f' {"--" if len(arg)>1 else "-"}{arg} {" ".join(kwargs[arg])}'
        elif type(kwargs[arg]) is bool:
            if kwargs[arg]:
                command += f' {"--" if len(arg)>1 else "-"}{arg}'
            else:
                continue
        else:
            command += f' {"--" if len(arg)>1 else "-"}{arg} {kwargs[arg]}'

    return command


def direct_call(command):
    '''
    Call a command directly.
    '''
    try:
        subprocess.check_call(command, shell=True)
    except subprocess.CalledProcessError as err:
        print(err)
        sys.exit(f"Error when running command '{command}'. Exiting.")
        pass
