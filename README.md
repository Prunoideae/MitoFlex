# MitoX

MitoX is a Python3 based toolkit mitochondrial genome assembling rewritten from [MitoZ](https://github.com/linzhi2013/MitoZ), with improved performance and result quality. And also for better extendability. It accepts both single-end and pair-end data, and follows a filter-assemble-annotate-visualize workflow to output results. Working mechanism is highly flexible and can be easily configured here.

# 1. System requirements

## 1.1 Platform

MitoX is developed under `Ubuntu 18.04.3 LTS on Windows Subsystem of Linux(WSL)`, and tested under `Ubuntu`

## 1.2 Storage

Installing MitoX requires abount 1~2GB of disk space. Assembling needs about 50% to 200% of raw data's size to create and store temp files each run. Temp files could be deleted after the run, though you may find something useful in it.

## 1.3 Memory

It takes about 5G to assemble the genome from a 10Gbp pair-end rawdata with thread number set to 80 (`--thread_number 80`). It takes much lower memory space in comparison to MitoZ, as MitoX uses the succinct de Brujin Graph (sDBG), a succinct representation of de Brujin Graph. Improving the data quality could reduce the memory usage.

## 1.4 CPU

MitoX makes use of a multi k-mer assemble strategy from [MEGAHIT](https://github.com/voutcn/megahit), thus requires more calculate power because multiple iteration of the graph is needed. Using MitoX with high iteration count (21-141, 8 generation) will takes more time than MitoZ in low thread number (8), but the time will be almost equivalent in a higher thread number (80).

## 1.5 GPU

MitoX does not explicitly requires GPU in the work, but a GPU will highly accelerate the process of sDBG building, about 3 times faster than normal.

# 2. Installation

## 2.1 From Docker Image

## 2.2 From git repository (Conda required)

# 3. Data requirement

# 4. Specifying parameters in configuration file

MitoX created a very flexible argument catching and processing mechanism, which is aimed to make it easier for further developing. [An example configuration file](example.config.py) is created under the main directory.

## 4.1 Configuration file structure

MitoX's configuration file is more like a independent python file than a traditional configuration file (*.cfg or *.json or something).

```python

# Most arguments are just the same like it's brother in commandline.
threads = 8

# The sub commands are specified by a recessive variable "command".
command = 'all'

# Also, you don't need to worry about the parameter order here.
# MitoX will just update everything correctly.

# Be aware: as the - character is not a valid character in the variables
# of Python, all the - are translated into _ in the scripts.
insert_size = 150

# Not specified parameters will use the default parameter, or parameters
# passed through commandline.

# MitoX first execute the configuration file as a independent python script,
# then collects all the variables under the global scope as parameters. So
# flexible parameter processing could be a possibility.
from os import path
fastq1 = path.abspath(
    path.join('~', 'data_folder', 'seq_paucumara_falcata', 'data.1.gz'))

# Also, MitoX will NOT check additional parameters passed through this way,
# only update and merge the parameters in config file with the commandline
# parameters. Temp variables will be passed but not processed.
foo = 'bar'

```

## 4.2 Generating configuration file

Besides creating a highly customized configuration file from the earth, you can also generates a configuration file from the commandline by specifying the option `--generate-config`, or `-g` in short. Then the MitoX will create a configuration file named `generated_config.py` under your current working directory.

Configuration generated in this way will have all the arguments needed by the specified subcommand (e.g. all or filter) in place, arguments passed latter will also be written into the config, overrides the default values predefined in scripts.

# 5. Modules provided by default MitoX

Most modules of MitoX are just the same as MitoZ. But some of the methods are rewritten or optimized.

## 5.1 all

## 5.2 filter

## 5.3 assemble

## 5.4 findmitoscaf

## 5.5 annotate

## 5.6 visualize

# 6. Extending the function of MitoX

I tried my best to make the code structure of MitoX as easy as possible to increase its extendability and readability for users to extend it if they find the tools used by MitoX are not good enough or the workflow could be optimized. Extending the function should not be an hard task as MitoZ now.

## 6.0 Calling the MitoX from other ways

Bash is not always the solution, in a certain circumstances an integrated call could be better because it allows deeper profiling and monitoring, and controlling.

Since the only goal of decorators is to expose the original methods of Python to the command line, calling the methods could be easy, but should be catious because a raw calling could bypass some of the fuse methods of MitoX, which could make the arguments parsed and processed not correctly, or even leads to an unknown end.

All MitoX function takes an Argument object as the only arguments, it's from the [parser.py](utility/parser.py), and it's very simple because actually it's just a wrapping of `**kwargs` to reduce coding and improve code redability.
Casting and passing the Argument instance is easy:

```python
from utility.parser import Arguments
args = Arguments({'foo':'bar', 'lorem':'ipsum'})
print(args.foo) # prints 'bar'
```

Once you obtain a Argument instance, you can directly call the function like this:

```python
import MitoX
MitoX.all(args)
```

But this is strongly NOT recommended, because it didn't follow the parser mechanism of MitoX, which makes the raw arguments to be passed to the function.

In here, a process function `process_arguments` is required from the [parser.py](utility/parser.py), it emulates how all the processer works and returns a processed Argument object, then you can call the subcommands in MitoX safely.

```python
from utility.parser import Arguments, process_arguments
import MitoX
# Assuming the groups and the handlers are all registered here.
args = Arguments({'foo':'bar'})
process_arguments(command='test', args=args)
MitoX.test(args)
```

## 6.1 Creating more subcommands

MitoX mainly use two decorators, `@parse_func` and `@arg_prop` from [parser.py](utility/parser.py), to profile and collect functions needed to be a subcommand of MitoX runtime, which means that  attaching the `@parse_func` decorator will expose the function to commandline, and using a `@arg_prop` decorator will add a argument to it.

```python
# Using a parse_func decorator will make this function 'visible' from command line.
# Also you can specify the parser groups to it.
@parse_func(func_help='example help', parents=[universal_parser])
# Using a arg_prop decorator will make this function tries to acquire an argument
# from command line, you can set many of the properties as you wish.
# By default you don't need to explicitly specify the type of variables, the MitoX
# will tried to use the type of default value at first, then turns to the str as
# it can accepts other types without losing data.
@arg_prop(dest='bar', help='fun args', default=1, type=float, required=True)
# An exception is the type bool, types or metavars will be disposed if the default
# value is a boolean, which makes the MitoX converts it to an action switch without
# considering other things. 'store_true' will be applied to variables with a False
# default value and vise versa.
@arg_prop(dest='switch', help='switchy', default=False)
# Another exception is the choices arguments, which will also disable specifying
# the type and metavar because all the type should be determined through the list.
@arg_prop(dest='list', help='of choices', choices=['foo', 'bar'], default='foo')
# The commandline underscore conversion is not applied here for the consistency
# of code.
@arg_prop(dest='under_score', help='haha')
# The function requires and only requires a object, Arguments, from the
# utility/parser.py, arguments will be passed as attributes of the object.
def foo(args):
    print(args.bar)
    print(args.under_score)

    # For further consideration, attributes are modifiable in the object, and
    # this will not change outside the scope because it's actually accessing
    # the class data. This could be powerful, but it needs you to take it
    # carefully.
    def bar(args):
        args.c = 1
    bar(args)
    print(args.c)
```

Calling the new added method from bash will be just like this and with output:

```bash
python3 MitoX.py foo --bar 2.0 --switch --under-score fun

2.0 # print(args.bar)
fun # print(under_score)
1   # print(args.c)
```

## 6.2 Creating parameter groups

Most methods shares a set of parameters, like thread numbers or input fastq file. Specifying the parameters repeatedly could be a problem, and validation will be difficult. So MitoX implements a parameter group processing mechanism to make this progress easier to be defined.
To create a argument group, you need to import and execute the method `register_group` from the module [parser.py](utility/parser.py):

```python
from utility.parser import register_group

# Argument groups support a callback function to reduce the code repeatence.
# You can modify and read the attributes straightforward and it will pass to the
# every downstream processors.(e.g. other handlers or the main function)
def handler(args):
    try:
        args.arglist = args.str_list.split(',')
    except Exception as i:
        print('Error occured when parsing the argument str-list!')
    # Returning value tells MitoX whether to run the process or not, returning
    # True means the processed arguments are valid, and vise versa.
        return False
    return True

foo_parser, foo_group = register_group('Test parser', [
    # Group arguments have almost the same options like the arg_prop function,
    # please refer to the arguments.py to check it out.
    {
        # The conversion rule applies here, because the arguments are quite
        # separated from the processing codes.
        # I much prefer - to _ because _ needs Shift + -, but - only needs one.
        'name':'str-list',
        'default':'1,2,3,4,5',
        'help':'input a set of numbers separated by comma(,)'
    }
], func=handler) # Specify the needed processor here.

# When multiple parser are specified, they follows an order of list parents to
# execute, this is important if there's a parser relies on other parser.
@parse_func(func_help='test func', parents=[universal_parser, foo_parser])
def bar(args):
    # Here we can use the argument directly created and processed by the handler.
    print(args.arglist)
```

Actually you can `register_group` anywhere as long as your `@parser_func` decorator can reach there, but I strongly recommends to write all the argument settings into the [arguments.py](arguments.py) to keep the code structure clean.