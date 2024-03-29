# MitoFlex

MitoFlex is a Python3 based toolkit designated for mitochondrial genome assembling, it's inspired from [MitoZ](https://github.com/linzhi2013/MitoZ), but with improved performance and result quality. And also it implemented a both easy and flexible mechanism to extend the program feature. It accepts both single-end and pair-end data, and follows an already set workflow to output results. Working mechanism is highly flexible and can be easily reconfigured here.

# 1 System requirements

## 1.1 Platform

MitoFlex is developed under `Ubuntu 18.04.3 LTS on Windows Subsystem of Linux(WSL2)`, compiled and tested under `CentOS release 7.3.1611`. Unix like system should work fine, but since some part of the program is compiled in `Ubuntu` or `CentOS`, MitoFlex may have risk to fail if running on other OS, like MacOS, Windows system is obviously not suitable to run MitoFlex, but `WSL` can, though there will be some performance loss.

## 1.2 File system

Installing MitoFlex requires about 1GB of space (Including dependency packages). Assembling needs about 50% to 200% of raw data's size to create and store temp files each run. Temp files could be deleted after the run, which will reduce the result to about 100MB size.

## 1.3 Memory

It takes about 5-30G to assemble the genome from a 5Gbps pair-end rawdata sample with thread number set to 80 (`--thread_number 80`). The memory consumption is highly varied, mainly depends on the fragmentation and the quality of rawdata, a dataset with more focused reads on mitogenome will absolutely takes much lesser memory since the contigs are limited. It takes much lower memory space in comparison to MitoZ, as MitoFlex uses the succinct de Brujin Graph (sDBG), a succinct representation of de Brujin Graph. Improving the data quality could reduce the memory usage. Also, more threads requires more memory in steps, since the data being processed is larger in the same time.

The average RAM consumption of MitoFlex is usually at 5G or even lower at 8 threads or 80 threads, but there are a few points that will consume much resource, the first one is the graph construction of megahit, which takes from 5G to 20G, but assembly will then only takes 1-3G to be done, the second one is the nhmmer search part explicitly in findmitoscaf module, where it takes about 15G or higher for searching against all the profiles. The steps are not parallelly processed, so 20GB or so is enough to deal with any type of data, or even less if the sample is small.

So, a machine with over 32 GB spare RAM is recommended, giving more could be more robust to deal with samples which are quite messy. As the workflow develops even further, some part of the module may requires even more or less RAM.

## 1.4 CPU

MitoFlex uses [megahit](https://github.com/voutcn/megahit) as assembler, thus requires more calculate power because multiple iteration of the graph is needed. The speed of assembly is mainly depends on how fragmentized the input reads are.

A recent modification enables MitoFlex to shorten the process time by over 75%, and improved the accuracy and sensitivity of the assembly. Using a dataset of 5Gbps ([Run browser](https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR1946581)), MitoFlex took 15min to finish the whole pipeline with 8 threads of Intel 9700KF and 20GB RAM usage on my laptop for developing.

## 1.5 GPU

MitoFlex does not explicitly requires GPU in the work, but a GPU will accelerate the process of sDBG building. The server I'm testing MitoFlex on has no GPU, so I can't tell much at this part.

# 2 Installation

## 2.1 From Docker Image (No recent implementation)

## 2.2 From git repository (Conda required)

For certain conditions, like if you don't have a sudo permission or root command, you can deploy MitoFlex from git directly.

To download MitoFlex from Github, simply type:

``` bash
git clone --depth=1 https://github.com/Prunoideae/MitoFlex
```

And git will pull the MitoFlex into your current directory, downloading the repo as zip and extract it to installation folder is also ok.

### 2.2.1 Installing Conda

Conda is required to create the environment MitoFlex needed instead, or you can install all the packages manually. Both [Anaconda](https://anaconda.org/anaconda/python) and [Miniconda](https://conda.io/miniconda.html) could be useful, but Miniconda is recommended if you don't need a big environment.

### 2.2.2 Setting up channels

These three channels should be the universal solution for any conda package installing.

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

As the network situation varies, mirror channels, or other alternative channels could be setup instead of the channels above for a better querying speed.

### 2.2.3 Creating environment for MitoFlex

Though you actually can setup the environment from the ground, using conda for creating virtual environment is recommended, as there may be other tools requiring different environment, and even requirements that are conflict to MitoFlex's, is actually a better choice.

MitoFlex requires a bunch of packages to run, you can install them in one line like:

```bash
conda create -n {environment name here} --file requirements.txt
```

All directly required packages are listed, if you don't want to use conda to install them:

```text
numpy pandas ete3 biopython psutil
megahit blast infernal circos hmmer wise2 bwa samtools infernal
```

You will have to solve the dependencies of required packages if not using conda. The upper 4 are python modules, and the lower ones are programs.

### 2.2.4 Setting up NCBI taxanomy database

Running `ncbi.py` from command line automatically updates the local database from NCBI taxanomy database.

```bash
./ncbi.py
```

Updating database from network is not always stable. `ncbi.py` will fall back to the local `taxdump.tar.gz` if any error occurs in the process.

### 2.2.5 Run MitoFlex

To test if MitoFlex is installed correctly, type :

```bash
./MitoFlex.py load_modules
```

If all modules are appeared to be correctly loaded, it indicates that you have all the required python modules installed, to have a test run on MitoFlex, please extract the files to somewhere from the `test.tar.gz` (in this [repository](https://github.com/Prunoideae/Samples)), and run :

```bash
/path/to/MitoFlex.py all workname test --use-list --fastq1 /path/to/fastq1 --fastq2 /path/to/fastq2
```

or you can use the `test_config.py` config file in the directory to test if the things are running correctly, just run:

```bash
/path/to/MitoFlex.py --config /path/to/config
```

this config assumes that the sample fastq is under the same directory as config does.

The test sample is done in 3min on my computer (Intel i7-9700KF, WSL2 Ubuntu) with 8 threads. 1-2 GB of spare RAM is required to run the sample.

Result will be the similar as that in `test.tar.gz` if there's no error in your installation, sequence should be cicular, PCGs and rRNAs should all be founded, otherwise you need to check the whole progress.

Exporting the directory to `PATH` environment variable is recommended for calling it more easily.

```bash
echo 'export $PATH="/path/to/installation/directory:$PATH"' >> '/path/to/rc'
```

Where rc stands for the .*rc file used by your terminal to load a certain set of commands at logging. Like `.bashrc` or `.zshrc`.

For a more detailed explanation of MitoFlex's arguments, please check out the program's help :

```bash
MitoFlex.py [module] <-h or --help>
```

This helps you better understand how MitoFlex will work, and so you can tune MitoFlex to meet your need in mitogenome assembly.

# 3 Data requirement

MitoFlex depends on the quality more than the size of data, it will not throw any error if your input fastq file is too small or something, but the result may be of low quality or nothing if the raw dataset is too small or unqualified to assemble a genome.

# 4 Specifying parameters in configuration file

MitoFlex uses a very flexible argument catching and processing mechanism, which is aimed to make it easier for further developing. [An example configuration file](example.config.py) is created under the main directory.

## 4.1 Configuration file structure

MitoFlex's configuration file is more like a independent python file than a traditional configuration file (\*.cfg or \*.json or something).

```python

# Most arguments are just the same like it's brother in commandline.
threads = 8

# The sub commands are specified by a recessive variable "command".
command = 'all'

# Also, you don't need to worry about the parameter order here.
# MitoFlex will just update everything correctly.

# Be aware: as the - character is not a valid character in the variables
# of Python, all the - are translated into _ in the scripts.
insert_size = 150

# Not specified parameters will use the default parameter, or parameters
# passed through commandline.

# MitoFlex first execute the configuration file as a independent python script,
# then collects all the variables under the global scope as parameters. So
# flexible parameter processing could be a possibility.
from os import path
fastq1 = path.abspath(
    path.join('~', 'data_folder', 'seq_paucumara_falcata', 'data.1.gz'))

# Also, MitoFlex will NOT check additional parameters passed through this way,
# only update and merge the parameters in config file with the commandline
# parameters. Temp variables will be passed but not processed.
foo = 'bar'

```

## 4.2 Generating configuration file

Besides creating a highly customized configuration file from the earth, you can also generates a configuration file from the commandline by specifying the option `--generate-config`, or `-g` in short. Then the MitoFlex will create a configuration file named `generated_config.py` under your current working directory.

Configuration generated in this way will have all the arguments needed by the specified subcommand (e.g. all or filter) in place, arguments passed latter will also be written into the config, overrides the default values predefined in scripts.

## 4.3 Specifying special parameters

Some parameters are freezed and specified to fit into the mitogenome analysis, they are either not performing good in most circumstances, or even will make the output result worse. But they may have some usage under special circumstances, so these parameters can be set in the [configuration](configurations.py) file.

# 5 Modules provided by MitoFlex

Most modules of MitoFlex are just the same as MitoZ. But some of the methods are rewritten and optimized.

## 5.1 all

Run the whole workflow, including methods listed below. Some part of module can be disabled to meet a certain usage, `all2` command is removed, and replaced with `--disable-filter` option, check the help for a more detailed usage.

## 5.2 filter

Filter out fastq sequences of low quality, binary is written in Rust to ensure speed and data safety. The method will not output compressed clean data by default, and most workflow is designed to directly process with plain data format, clean data will be deleted after the `all` command is done if `--keep-temp` option is not set.

## 5.3 assemble

Assemble the fastq file to output contigs. This method use megahit for faster and better results, but since megahit itself implemented a multi-kmer strategy to assemble data, it might take long (average 20min for a 5Gbps dataset, ranging from 1min to 1.5h if filtered correctly with 8 threads, will be much more faster if more processors specified) to assemble. Reducing kmer steps, or disabling local assembly could shorten the time, but it's not recommended since it also may output more fragmentized contigs. Also, assemble process depends much more on data quality than the size of dataset, because it will take much more resources to process more contigs in each iteration, if final sequence itself is fragmentized.

To eliminate this problem and make the assembler stay focus on the mitogenome sequences, MitoFlex tweaked how the megahit's original pipeline works, inserted a filter method to filter out sequence of not enough depth - which is considered to be nuclear genome or contamination, since mitogenome should possess a much higher copy number in the reads. This strategy can be configured via `--depth-list` option. If you need to disable the option to find out some sequence, like NUMTs, please check out `configurations.py`.

## 5.4 findmitoscaf

To pick out candidate sequences which likely to be mitochondrial sequences. This process firstly drop sequences which is obviously too low for searching, then will search the data with nhmmer and tblastn. The two step of process, nhmmer and tblastn, depends highly on the parallel processing power of the environment, in our testing environment (80 threads), this takes 5-10min for any dataset to finish.

## 5.5 annotate

To annotate sequences using tblastn and infernal.

## 5.6 visualize

To generate PNG and SVG file representing current mitogenome in a more direct way.

# 6 Frequently occured problems

## 6.1 "Killed" or showing not enough memory, or something like that

Your machine is running out of memory, please add more. This is not a code bug or something, it's necessary for MitoFlex to use some memory, to assemble from GBs of reads, or picking up candidate from millions of contigs.

## 6.2 Mitogenome sequences not found in `findmitoscaf`, or not outputting enough sequences in assembly

There are several possible problems to make this happen. Check if your run is :

1. of low bps? Like only 1.5Gbps of the data, if so, please specify a more tolerant depth list, less reads decrease the overall depth of sequences, your result may be filtered out during the iteration.
2. of low quality? It could be happened sometimes, if your raw data is of too many other sequences, like contaminated or just because it's poor, then it will happens as above.
3. specified a too strict depth list? In most test cases, target sequence depths varies from 300 to 1000, in dedicated result it can go even up to 10000, result sequences will just be dropped out if your sequences is of not enough depth.
4. analysing species of small database? `findmitoscaf` module depends on database much, if the species you want to assemble is having too few data, it could make MitoFlex fail to detect mitogenome sequences, if so, please extend the database by adding more profile.

## 6.3 "Iteration broke at kmer = x, since no no valid contig in kmer = y is done!"

This is because no sequences is output during the last iteration, no any of the contig will output if the iteration continues, so MitoFlex terminate the assembly. It's similar to problem `6.2`.

## 6.4 Input file 1 and 2 have different sizes! This could cause loss on rawdata, or even crash the program

Your input PE reads are not of the same size, it indicates a mismatch between two file, and it will likely stop the programs like `bwa` from starting, halting the pipeline, or have some graph element missing. In most of time maybe this will work fine, but there maybe times that later program refuse to startup.

## 6.5 Module XXX is not found! / Cannot import helper module XXX

You have an invalid installation, or environment setup, please check if your installation is good, and you have switched to a proper conda environment.

## 6.6 Cannot validate folder / Folder X already exists

MitoFlex failed to create a folder, or a previous folder exists so MitoFlex can't make sure if the folder is removable or not. Please check your permission or try to investigate the files in the folder and remove it if it's safe to remove.

## 6.7 I'm using MitoFlex on flatworms, but the closely-related species is in Echinodermata

MitoFlex doesn't know what species you are processing, the only thing it has just a local protein database, which contains many protein records of different species.

Since then, such a species-judging method is not driven by something that so accurate, MitoFlex can only check out the numbers of PCGs that certain species contains, if a species contributed the most PCGs during the annotation, then it's the closet-related progress, and if two or more species contributed the same amount, only a random *one* of them can be selected.

So, the accuracy is neither guaranteed nor possible. If you want to use MitoFlex to aid your classification of species, I would suggest you to check out the `wise.csv` under the `workname/temp/annotation` folder, it contains a full list of PCGs and corresponding species that are thought to be the mitogenome's hit, you may find something more in this.

## 6.8 Running Circos error

Your genome sequences is oversized, even after some deduplication it still have over 25000 points to draw, which is too many for drawing a good and clean plot.

Such a large genome should not occur commonly, if you can sure about this genome is what you want, please change your `max_points_per_track` in Circos settings to a larger number to make Circos actually draw this.

Such a problem only affects visualization, no effect on previous methods.

## 6.9 Processing species not included in database

If there's enough support data, you can refer to chapter 7 to add data to build a taxonomy database for the species. Any pull requests to add new taxonomy clade is gracefully welcomed.

Or, you can just use a clade close to the species you have. In this way, filtering, assembling and mitogenome scaffolding merging and finding will not be influenced much, while annotation could have problem since it depends heavily on the local database. External usage of annotation tools, like `MITOS` or others is then recommended in this situation.

# 7 Adding new profile data to MitoFlex

MitoFlex has already integrated protein and nucleic acid data into the profile, but it can't cover all the species for sure. So it's necessary to add data of other taxonomy classes if current MitoFlex doesn't have it for better assemble and annotation performance.

## 7.1 Building nhmmer profile

For official documentation refers to [here](http://www.csb.yale.edu/userguides/seq/hmmer/docs/node19.html).

Building HMM profile needs the access of `hmmbuild`, which is included in the HMMER package, the command requires a Multiple Sequence Alignment (MSA) file, which is obtained from aligning sequences you want to build profile with by serveral alignment program like MAFFT or ClustalW. The `hmmbuild` in installed version supports most file format, like FASTA, Stockholm or ClustalW.

## 7.2 Adding profile data

MitoFlex has its internal profile for basic mitogenome assembly and annotation, but it comes to be inaccurate if a specific speciemen is required. You can implement your own set of profile of the specific speciemen you want if feeling like MitoFlex is not giving good results.

### 7.2.1 Adding or modifying clade protein database

MitoFlex uses a given set of protein sequences to do tblastn, for picking up most possible sequences and for annotating the genes, if you found the species you requested is not quite covered in the database (Like Rhabditophora in Playthelminthes), you can of course add your sequences into the profile. The set of sequences can be found in `profile/MT_database/{clade}.fa`, written in FASTA file format, the id of the sequence must be in `gi_NC_{record id}_{gene}_{genus}_{species}_{length}_aa`, only gene, genus and species were taken into recognition of sequences and clades, but please keep underscores in place for the program to detect and parse the information.

### 7.2.2 Adding a new clade

To add a new clade, these three files should be noticed: `{clade}.hmm` and `required_cds.json` in `profile/CDS_HMM` and `{clade}.fa` in `profile/MT_database`. The hmm file is used for nhmmer to search out possible sequences, and the fasta file is for tblastn to mark the potential sequence with genes, and the `required_cds.json` is for telling MitoFlex what gene should be taken into count, because some gene is rarely reported in several species, so it would be better for users to tell how MitoFlex will judge the gene is missing or not. MitoFlex requires all these three to be set properly in order to be functional.

### 7.2.3 Adding Covariance Models for tRNA search

Please put your cm file into the [tRNA_CM](profile/tRNA_CM) folder, MitoFlex will automatically use files under this directory for tRNA searching.

# 8 Things that will effect MitoFlex's overall performance

MitoFlex itself isn't doing magic, though it's fast and (more) reliable, bad result could be outputted if you don't even give it a good input, but at least it will be the most reliable result of all contigs that megahit can assemble, and it will be better than MitoZ.

There are mainly two factors influencing result quality : 1. The quality of rawdata itself. 2. The size of genome profile that MitoFlex currently have in the profile folder.

## 8.1 Data quality

MitoFlex doesn't depends on the size of data too much, since only a portion of total reads is filtered and extracted, resulted to be the clean data, the data actually used in the workflow is small, no more than 5Gbps(can be adjusted if needed). The speed of assembly is then depends on the data quality, which is mainly representing how filtered reads are concentrated on your final sequence, this is quite obvious, if the reads are aparted and have no linkage, the assembler will have to retain them at each iteration because you can't actually concat them into a single sequence, which leads to a slower iteration, and a bigger contig file. So, the findmitoscaf will have to pick out mitogenome sequences from a much larger contig file, where it could also be quite slow.

The other part of data quality, is how much mitogenome is covered in the raw data, this could be varied a lot, if the input data is small, or just because the filter module extracted too much nuclear genomic sequences. Since the mitogenome is really small and possess a high depth number over the nuclear genome, the latter is rarely happened. But if you ensures that your dataset is of enough size and quality, you can increase the truncation threshold by specifing `--trimming <INT>` to X Gbps you want.

## 8.2 Genome profile

MitoFlex uses `nhmmer` and `tblastn` to search for mitogenome candidates, where `nhmmer` is used to identify how a region on a sequence is related to some PCG, which is quite remote but accurate enough. This process is of high tolerance and profile can be used across phylums. BUT, the accuracy, and the average length of alignment, will be then severely reduced, where it may be difficult for MitoFlex to identify the fragments of contigs, or to determine where the sequence is belong to.

The `tblastn` is used for a more closely related homology search, and used to determine and filter out the sequences not belonging to the taxanomy clade expected, which is categorized as contamination. But if none of the gene related to the clade annotated by searching through `tblastn`, the program will remove it incorrectly. Also it will affect the `annotation` module, so there may be a mismatch for the `findmitoscaf` and the `annotation` results, where mainly indicates the protein database currently have is of not enough records, and an additional search using other methods of annotation, for example a `MITOS` web server is recommended.

Further adjustment on the process may be scheduled, but not now.

# 9 Extending the function of MitoFlex

Although MitoFlex has already implemented a full workflow to filter, assemble and annotate the mitogenome, and all of this can be done in one-click, it also support to modify some behaviour if you want to do. This section is for extending MitoFlex for your own usage, though most users may not have the need to do this.

MitoFlex is designed for extendability and readability, to make users to extend it if they find the tools used by MitoFlex are not good enough or the workflow could be even optimized. Extending the function should not be an hard task as MitoZ.

MitoFlex is written in Python 3.6, so modifying the original workflow of MitoFlex requires a basic knowledge of the Python programming language.

## 9.0 Calling the MitoFlex from other ways

Bash is not always the solution, in a certain circumstances an call from python inside could be better because it allows deeper profiling and monitoring, and controlling. For example deploying and integrating the MitoFlex with environments like `Jupyter` or web servers like `Django` or `Flask`.

Since the only goal of decorators is to expose the original methods of Python to the command line, calling the methods could be easy, but should be catious because a raw calling could bypass some of the fuse methods of MitoFlex, which could make the arguments parsed and processed not correctly, or even leads to an unknown end.

All MitoFlex function takes an Argument object as the only arguments, it's from the [parser.py](utility/parser.py), working like a context argumentin certain CLI builder like `invoke` or others, it's very simple because actually it's just a wrapper of `**kwargs` to reduce coding and improve code redability.
Casting and passing the Argument instance is easy:

```python
from utility.parser import Arguments
args = Arguments({'foo':'bar', 'lorem':'ipsum'})
print(args.foo) # prints 'bar'
```

Once you obtain an Argument instance, you can directly call the function like this:

```python
import MitoFlex
MitoFlex.all(args)
```

But this is strongly NOT recommended, because it didn't follow the parser mechanism of MitoFlex, which makes the raw arguments to be passed to the function.

In here, a process function `process_arguments` is required from the [parser.py](utility/parser.py), it emulates how all the processer works and returns a processed Argument object, then you can call the subcommands in MitoFlex safely.

```python
from utility.parser import Arguments, process_arguments
import MitoFlex
# Assuming the groups and the handlers are all registered here.
args = Arguments({'foo':'bar'})
process_arguments(command='test', args=args)
MitoFlex.test(args)
```

## 9.1 Creating more subcommands

MitoFlex mainly use two decorators, `@parse_func` and `@arg_prop` from [parser.py](utility/parser.py), to profile and collect functions needed to be a subcommand of MitoFlex runtime, which means that  attaching the `@parse_func` decorator will expose the function to commandline, and using a `@arg_prop` decorator will add a argument to it.

```python
# Using a parse_func decorator will make this function 'visible' from command line.
# Also you can specify the parser groups to it.
@parse_func(func_help='example help', parents=[universal_parser])
# Using a arg_prop decorator will make this function tries to acquire an argument
# from command line, you can set many of the properties as you wish.
# By default you don't need to explicitly specify the type of variables, the MitoFlex
# will tried to use the type of default value at first, then turns to the str as
# it can accepts other types without losing data.
@arg_prop(dest='bar', help='fun args', default=1, type=float, required=True)
# An exception is the type bool, types or metavars will be disposed if the default
# value is a boolean, which makes the MitoFlex converts it to an action switch without
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
    print(args.switch)

    # For further consideration, attributes are modifiable in the object, and
    # this will not be dropped outside the scope because it's actually accessing
    # the class data. This could be powerful, but it needs you to take it
    # carefully.
    def bar(args):
        args.c = 1
    bar(args)
    print(args.c)
```

Calling the new added method from bash will be just like this and with output:

```bash
python3 MitoFlex.py foo --bar 2.0 --switch --under-score fun

2.0     # print(args.bar)
fun     # print(args.under_score)
True    # print(args.switch)
1       # print(args.c)
```

## 9.2 Creating parameter groups

Most methods shares a set of parameters, like thread numbers or input fastq file. Specifying the parameters repeatedly could be a problem, and validation will be difficult. So MitoFlex implements a parameter group processing mechanism to make this progress easier to be defined.
To create a argument group, you need to import and execute the method `register_group` from the module [parser.py](utility/parser.py):

```python
from utility.parser import register_group

# Argument groups support a callback function to reduce the code repeatence.
# You can modify and read the attributes straightforward and it will pass to the
# every downstream processors.(e.g. other handlers or the main function)
def handler(args):
    try:
        args.arglist = args.str_list.split(',')
    except Exception:
        print('Error occured when parsing the argument str-list!')
    # Returning value tells MitoFlex whether to run the process or not, returning
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

# 10 Reusing my code

I'm very glad to see that my code is used in other fields, even in non-bioinformatic way. But sadly this toolkit is done in a hurry, and lacks detailed utilization of design patterns, so some part of it may look ugly, and even confusing to understand. This toolkit have incubated many other repositories that I've found something worth to deal with, but they all require some time before being developed, since I have to finish my courses as a undergraduate in SZU.

The [utility](utility/) folder contains most helper classes and methods used in the program development. If you want to directly reuse the code I wrote in related research field(mitogenome analyzing, for example), please cite my paper if you will publish one, if you want to implement a similar workflow yourself, please cite the MitoZ's paper since this program is inspired from the former toolkit.

The argument [parser](utility/parser.py) of this software is quite useful, but sadly I strongly not recommend you to use it directly in the program, since it was a temporarily made argparse wrapper just in three days, much code here is neither clean, nor easy to use, though already suitable for MitoFlex's current need and hard to rewrite one in short time. If you want to implement a workflow into a similar framework without having too much coding(like directly facing the argparse module), better integration of other caller besides bash or command prompt, for example a `Jupyter` web notebook or a `Django` or `Flask` server, and more flexible, modularized, clean code and code structure, please refer to the `ArgWraps`, a [repo](https://github.com/Prunoideae/ArgWraps) pinned on my GitHub page.

For the analyzer of Washington University Secondary Structure (WUSS), any rewriting is welcomed, since the output format of Infernal is quite messy, the parser of this can only accept annotations in one line (though it can parse any structure if sequence and fold string are given directly), and I have no idea on how to parse the file since it seems to have no actual format for me to write the parser. Any help on this is gratefully welcomed.
