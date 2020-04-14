"""
filter.py
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
from os import path
try:
    sys.path.insert(0, os.path.abspath(os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..")))
    from utility.helper import shell_call
    from utility import logger
except ImportError as err:
    sys.exit(f"Unable to import helper module {err.name}, is the installation of MitoFlex valid?")

filter_dir = os.path.dirname(os.path.abspath(__file__))


def filter_se(fqiabs=None, fqoabs=None, Ns=10, quality=55, limit=0.2, start=None, end=None, trim=0, trunc=False):
    fsin = path.getsize(fqiabs)
    logger.log(level=1, info='Start filtering single-end rawdata.')
    logger.log(level=0, info=f'Input file has {fsin} bytes.')
    logger.log(level=1,
               info=f'Using argument : Ns={Ns}, quality={quality}, limit={limit}, start={start}, end={end}, trimming={trim}, trunc={trunc}')
    try:
        shell_call(path.join(filter_dir, 'filter_v2'), cleanq1=f'"{fqoabs}"', fastq1=f'"{fqiabs}"',
                   n=Ns, q=quality, l=limit, s=start, e=end, t=trim, truncate_only=trunc)
    except Exception as identifier:
        logger.log(
            level=4, info=f'Error occured when running filter, cause : {identifier}')
        logger.log(level=1, info=f'Input file : {fqiabs}')
        logger.log(level=1, info=f'Output file : {fqoabs}')

        sys.exit("Error occured when running filter!")

    fsot = path.getsize(fqoabs)
    logger.log(level=0, info=f'Output file has {fsot} bytes.')
    logger.log(level=0,
               info=f'Filtered {fsin - fsot} bytes, ratio {fsot/fsin}.')

    return fqoabs


def filter_pe(fq1=None, fq2=None, o1=None, o2=None,
              dedup=False, start=None, end=None,
              n=10, q=55, l=0.2, trim=0, trunc=False):
    fsin1, fsin2 = path.getsize(fq1), path.getsize(fq2)
    logger.log(level=1, info='Start filtering pair-end rawdata.')
    logger.log(
        level=0, info=f'Input file 1 has {fsin1} bytes, 2 has {fsin2} bytes.')
    if fsin1 != fsin2:
        logger.log(
            level=3, info=f'Input file 1 and 2 have different sizes! This could cause loss on rawdata, or even crash the program.')
    logger.log(
        level=1, info=f'Using argument : Ns={n}, quality={q}, start={start}, end={end},limit={l}, trimming={trim}')
    try:
        shell_call(path.join(filter_dir, 'filter_v2'),
                   _1=f'"{fq1}"', _2=f'"{fq2}"', _3=f'"{o1}"', _4=f'"{o2}"', d=dedup, s=start,
                   e=end, n=n, q=q, l=l, t=trim, truncate_only=trunc)
    except Exception as identifier:
        logger.log(
            level=4, info=f'Error occured when running filter, cause : {identifier}')
        logger.log(level=1, info=f'Input file : {fq1} , {fq2}')
        logger.log(level=1, info=f'Output file : {o1} , {o2}')
        sys.exit("Error occured when running filter!")

    fsot1 = path.getsize(o1)
    logger.log(level=0, info=f'Output file has {fsot1} bytes.')
    logger.log(level=0,
               info=f'Filtered {fsin1 - fsot1} bytes, ratio {fsot1/fsin1}.')
    return o1, o2
