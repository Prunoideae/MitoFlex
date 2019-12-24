"""
profiler.py
=========

Copyright (c) 2019-2020 Li Junyu <2018301050@szu.edu.cn>.

This file is part of MitoX.

MitoX is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MitoX is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MitoX.  If not, see <http://www.gnu.org/licenses/>.

"""

import inspect
import time
from os import path
import os
import sys
import json
from threading import Thread

try:
    import psutil
    import shutil
except:
    print('Error occured when importing psutil, profiler will be disabled.')

# Using matplotlib to visualize the data, not necessary to the
# profiling.
plot_enabled = True
try:
    import matplotlib
except ImportError as ident:
    plot_enabled = False


# Inner side profiler.
class Profiler():
    '''
    A wrapper class to record performance usage of certain phase.
    '''

    def __init__(self, tickrate=1):
        self._runner = Thread(target=__loop, args=[self])
        self._state = 0
        self.tickframes = []
        self.tickrate = tickrate

    def __loop(self):
        base_time = time.time()
        while self._state != 2:
            self.tickframes.append(
                {
                    'time': time.time() - base_time,
                    'cpus': psutil.cpu_percent(percpu=True),
                    'memory': psutil.virtual_memory(),
                    'disks': psutil.disk_io_counters(perdisk=True)
                }
            )
            time.sleep(self.tickrate)

    def start(self):
        if self._state != 0:
            raise RuntimeError('The profiler is already started or stopped.')
        self._state = 1
        self._runner.start()

    def stop(self):
        self._state = 2
        pass


def profiling(func):
    '''
    Quest a simple profiler for fetching data of running every method.
    This module is created to collect statistical data for further
    analysis and optimizing.

    Needs to pass a additional argument 'profiling_path' to the
    injected method to activates the profiling, or the decorator
    will never understand where to write reports.

    Passing None to the profiling_path will also cause the profiler
    not to be running, but be noticed, this profiler will NOT check
    if the profiling path is valid. Handling should be done as the
    parser regulator be called.

    Also, the profile decorator will shield the decorated function
    from exposing, as the profiler is actually a wrapper sugar.

    This decorator should not be applied to any exposed function, as
    it will cause trouble for obvious reason.'''

    def wrapper(*args, **kwargs):

        profiling_path = None
        profiler = None
        if 'profiling_path' in kwargs:
            profiling_path = kwargs['profiling_path']
            kwargs.pop('profiling_path')
            profiler = Profiler(
                tickrate=1 if 'tickrate' not in kwargs else kwargs['tickrate'])
            profiler.start()

        if profiling_path is not None:
            profiling_path = path.abspath(profiling_path)
            try:
                os.makedirs(profiling_path, exist_ok=True)
            except Exception as identifier:
                sys.exit('Error occured when validating directories.')

        func(*args, **kwargs)

        if profiling_path is not None:
            profiler.stop()
            with open(path.join(profiling_path, func.__name__ + '.stat'), 'w') as f:
                json.dump(obj=profiler.tickframes, fp=f)

    return wrapper


def collect_results(profile_dir=None, plotting=True):
    if plotting and not plot_enabled:
        print("Matplot not enabled, only text summary will be enabled.")
    plotting = plotting and plot_enabled
