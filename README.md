# MitoX

MitoX is a Python3 toolkit rewritten from [MitoZ](https://github.com/linzhi2013/MitoZ), a mitochondrial genome assembling toolkit to imporve the performance and result quality. And also for imporving the extendability. It accepts both single-end and pair-end data, and follows a filter-assemble-annotate-visualize workflow to output results.

# 1. System requirements

## 1.1 Platform

MitoX is developed under `Ubuntu 18.04.3 LTS on Windows Subsystem of Linux(WSL)`, and tested under `Ubuntu`

## 1.2 Storage

Installing MitoX requires abount 1~2GB of disk space. Assembling needs about 50% of raw data's size to create and store temp files each run. Temp files could be deleted after the run, though you may find something usefule in it.

## 1.3 Memory

It takes about 5G to assemble the genome from a 10Gbp pair-end rawdata with thread number set to 80 (`--thread_number 80`). It takes much lower memory space in comparison to MitoZ, as MitoX uses the succinct de Brujin Graph (sDBG), a succinct representation of de Brujin Graph. Improving the data quality could reduce the memory usage.

## 1.4 CPU

MitoX makes use of a multi k-mer assemble strategy from [MEGAHIT](https://github.com/voutcn/megahit), thus requires more calculate power because multiple iteration of the graph is needed. Using MitoX with high iteration count (21-141, 8 generation) will takes more time than MitoZ in low thread number (8), but the time will be almost equivalent in a higher thread number (80).

## 1.5 GPU

MitoX does not explicitly requires GPU in the work, but a GPU will highly accelerate the process of sDBG building, about 3 times faster than normal. As the calculation could be done parallely

# 2. Installation

## 2.1 From Docker Image

## 2.2 From git repository

# 3. Data requirement

# 4. Specifying parameters in configuration file

MitoX created a very flexible argument catching and processing mechanism, which is aimed to make it easier for further developing. [A example configuration file](example.config.py) is created under the main directory.

## 4.1 Configuration file structure

MitoX's configuration file is more like a dependent python file to a traditional configuration file (*.cfg or *.json or something).

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
# only update and merge the parameters in config file into the commandline
# parameters, so use temp variable as you favor.
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

## 6.1 Creating more subcommands

MitoX mainly use two decorators, `parse_func` and `arg_prop` from [parser.py](utility/parser.py), to profile and collect functions needed to be a subcommand of MitoX runtime, which means that  attaching the `parse_func` decorator will expose the function to commandline, and using a `arg_prop` decorator will add a argument to it.

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
