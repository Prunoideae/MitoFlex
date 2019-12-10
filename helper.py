"""
helper.py
=========

Copyright (c) 2019-2020 Henry Lee <2018301050@szu.edu.cn>.

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

import argparse
import inspect
import subprocess
import sys
import os
from functools import wraps, partial
from importlib import util as import_util

collected_args = {}

# Argument processing


def register_group(group_name, argument_list):
    '''
    Create an argument group with group_name and dict argument_list. 
    Returns parser and group of this.
    '''
    parser = argparse.ArgumentParser(add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_argument_group(group_name)

    for i in argument_list:
        # Determine which type will be used in the argument,
        # default to str, but will be changed if type or default are set.
        real_type = str
        if 'type' in i:
            real_type = i['type']
        elif 'default' in i:
            real_type = type(i['default'])

        # Use a tricky way to generate metavars for arguments.
        real_meta = i['meta'] if 'meta' in i else real_type.__name__

        real_help = i['help'] if 'help' in i else ''
        real_help += f" Default {i['default']}." if 'default' in i and i['default'] != '' and real_help != '' else ''

        if 'action' in i:
            group.add_argument(
                f'--{i["name"]}',
                required=i['required'] if 'required' in i else False,
                action=i['action'],
                default=i['default'] if 'default' in i else None,
                help=real_help)
        elif 'choices' in i:
            group.add_argument(
                f'--{i["name"]}',
                required=i['required'] if 'required' in i else False,
                choices=i['choices'],
                default=i['default'] if 'default' in i else None,
                help=real_help,
                type=real_type)
        else:
            group.add_argument(
                f'--{i["name"]}',
                metavar=f"<{real_meta.upper()}>",
                required=i['required'] if 'required' in i else False,
                default=i['default'] if 'default' in i else None,
                help=real_help,
                type=real_type)

    return parser, group


def parse_func(func=None, *, func_help='', parents=[]):
    '''
    Mark the decorated function as a valid 'argument acceptable' 
    function. Decorated function will be analysed once it's created, 
    and it will further be examined and called.\n
    Values' metavars and types are determined automatically unless it
    was modified before the parse_func decorator, but be careful, every 
    bool argument will NOT be treated as a normal argument but a switch, 
    thus you can't specify a metarvar or assign a type to it. Though modify 
    the argument attributes after the decoration is acceptable, it's 
    not recommended for code structure.\n
    Arguments introduction:\n
    func : a preserved variant for further call, DO NOT MODIFY IT!\n
    func_help : give the function a help string, which will be shown
    in --help or parser.print_help()\n
    parents : specify the parent argument parsers for this func.
    '''
    if func is None:
        return partial(parse_func, func_help=func_help, parents=parents)

    global collected_args

    func_name = func.__name__
    if func_name not in collected_args:
        collected_args[func_name] = {
            'args': {}
        }
    collected_args[func_name]['valid'] = True
    collected_args[func_name]['help'] = func_help
    collected_args[func_name]['parents'] = parents
    collected_args[func_name]['func'] = func

    shadowed = []
    if parents is not None:
        for parent in parents:
            for group in parent._action_groups:
                shadowed += [x.dest for x in group._group_actions if x.dest not in shadowed]

    spec = inspect.getargspec(func)
    func_args = collected_args[func_name]['args']
    if spec.defaults is not None:
        for arg, default in zip(spec.args, spec.defaults):
            if arg not in func_args:
                func_args[arg] = {
                    'type': type(default) if default is not None else str,
                    'default': default if default is not None else '',
                    'help': '',
                    'required': False,
                    'choices': None
                }

    for arg in func_args.copy():
        if arg in shadowed:
            func_args.pop(arg)

    return func


def arg_prop(func=None, *, dest=None, arg_type=None, help=None, required=None, choices=None, meta=None):
    '''
    Specify a destination var, then modify its attributes.\n
    Argument introduction:\n
    func : a preserved variant for further call, DO NOT MODIFY IT!\n
    dest : specify the destination argument, leaving this to blank or specify 
    a argument that's not existed will raise an error.\n
    arg_type : specify the argument type of the argument, like you may change 
    int to float, or something alike. A argument have a boolean default will not be changed.\n
    help : add a help string to argument. This will appear in --help or parser.print_help().\n
    required : mark this argument is required or not.\n
    choices : give the argument a certain choices, this is predefined and not immutable, neither 
    it has a metavar nor undecleard value can be entered in this argument.\n
    meta : specify the metavar of the argument, like str can be changed to file, for a more detailed 
    information, this will not influence how the argument works.\n
    '''
    if func is None:
        return partial(arg_prop, dest=dest, help=help, required=required, choices=choices)

    if dest is None:
        raise AttributeError('No arg specified to change properties!')

    global collected_args

    func_name = func.__name__
    if func_name not in collected_args:
        collected_args[func_name] = {
            'valid': False,
            'help': '',
            'args': {}
        }

    spec = inspect.getargspec(func)
    func_args = collected_args[func_name]['args']
    if dest not in spec.args:
        raise AttributeError(
            f"Argument {dest} not found in function {func_name}!")
    if dest not in func_args:
        default = spec.defaults[spec.args.index(dest)]
        func_args[dest] = {
            'type': arg_type if arg_type is not None and type(default) is not bool
            else (type(default) if default is not None else str),
            'default': default if default is not None else '',
            'help': help if help is not None else '',
            'required': required if required is not None else False,
            'choices': choices,
            'meta': meta
        }
    else:
        arg = func_args[dest]
        if arg_type is not None and arg['type'] is not bool:
            arg['type'] = arg_type
        if help is not None:
            arg['help'] = help
        if required is not None:
            arg['required'] = required
        if choices is not None:
            arg['choices'] = choices
        if meta is not None:
            arg['meta'] = meta
    return func


def freeze_arguments(prog, desc):
    '''
    Collects all the information needed, and create a universal parser for all sub parser.
    This will not lock the global argument data, so multiple parser may be created seperately.\n
    Argument introduction:\n
    prog : what is the program.\n
    desc : add a description for the parser.
    '''
    global collected_args
    main_parser = argparse.ArgumentParser(
        prog=prog, description=desc, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = main_parser.add_subparsers(dest='command')
    for f in collected_args:
        func = collected_args[f]
        if not func['valid']:
            continue
        parser_func = subparsers.add_parser(
            f, parents=func['parents'], help=func['help'])
        for arg in func['args']:
            arg_attrs = func['args'][arg]
            if arg_attrs['type'] is bool:
                parser_func.add_argument(
                    f'--{arg}',
                    default=arg_attrs['default'],
                    action='store_false' if arg_attrs['default'] else 'store_true',
                    help=arg_attrs['help'],
                    required=arg_attrs['required']
                )
            elif arg_attrs['choices'] is not None:
                parser_func.add_argument(
                    f'--{arg}',
                    default=arg_attrs['default'],
                    choices=arg_attrs['choices'],
                    help=arg_attrs['help'],
                    required=arg_attrs['required']
                )
            else:
                parser_func.add_argument(
                    f'--{arg}',
                    default=arg_attrs['default'],
                    help=arg_attrs['help'],
                    required=arg_attrs['required'],
                    metavar=f"<{arg_attrs['type'].__name__.upper() if arg_attrs['meta'] is None else arg_attrs['meta']}>",
                    type=arg_attrs['type']
                )
        func['parser'] = parser_func
    main_parser.add_argument("-c", "--config", type=str, metavar='<FILE>',
                             help='use preconfigurated file to run program')
    main_parser.add_argument("-g", "--generate_config", type=bool, action="store_true", default=False,
                            help=("if switched on, MitoX will not be run, but generate"
                                "a configuration file with arguments input instead."))
    return main_parser


def parse_then_call(expr):
    '''
    Analyze the parser information, then call a certain function with the name 
    registered with parse_func before.\n
    Argument introduction:\n
    expr : Accepts a parser created from freeze_argument. Other parsers could 
    be used, but unregistered function will make the program extremely unstable 
    and may lead to an unhappy end.\n
    '''
    args = expr.parse_args()
    parsed = vars(args)

    config = parsed['config']
    parsed.pop('config')

    if config is not None:
        try:
            spec = import_util.spec_from_file_location('', config)
            mod = import_util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            parsed.update({
                item: vars(mod)[item]
                for item in vars(mod) 
                if not item.startswith("__")
            })
        except Exception as identifier:
            print(identifier)
            sys.exit('Errors occured when importing configuration file. Exiting.')

    command = parsed['command']
    if command in collected_args:
        parsed.pop('command')
        func = collected_args[command]['func']
        func(**parsed)
    else:
        expr.print_help()

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
    '''
    command = ' '.join(args)
    for arg in kwargs:
        if type(kwargs[arg]) is not list:
            command += f' --{arg} {kwargs[arg]}'
        else:
            command += f' --{arg} {" ".join(kwargs[arg])}'
    direct_call(command)


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
