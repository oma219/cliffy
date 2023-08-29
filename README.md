# **Document Array Profiles**: document listing in compressed space [![DOI](https://zenodo.org/badge/535600494.svg)](https://zenodo.org/badge/latestdoi/535600494)

This code is based on the [pfp-thresholds](https://github.com/maxrossi91/pfp-thresholds) repository written by [Massimiliano Rossi](https://github.com/maxrossi91). 
This repository **builds a novel data-structure called *document array profiles* that allows users to list all the documents/classes that
a particular exact match occurs**, theoretically described by Massimiliano Rossi. This data-structure is built based on the prefix-free parse of the 
input text.

## Installation

For starting out, use the commands below to download the repository and build the executable. After running the make command below,
the `pfp_doc` executable will be found in the `build/` folder.

```sh
git clone git@github.com:oma219/docprofiles.git
cd docprofiles

mkdir build 
cd build && cmake ..
make install

export PFPDOC_BUILD_DIR=$(pwd)
```

## Getting started

The basic workflow with `pfp_doc` is to (1) build the data-structure over an input text and then (2) run queries
using the document array profiles.

### (1) Build the data-structure 

```sh
pfp_doc build  -f <input_list> -o <output_prefix>
```

The command above shows how to build the doument array profiles data-structure, it takes in a file-list of multiple genomes
and then generates an output file/index using the output prefix. In the file-list, you can specify a list of 
genomes and then specify which document/class each genome belongs in.

**Example of file-list file:**
```sh
/path/to/ecoli_1.fna 1
/path/to/ecoli_2.fna 1
/path/to/salmonella_1.fna 2
/path/to/human_1.fna 3
/path/to/human_2.fna 3
```

### (2) Run queries

```sh
pfp_doc run -r <output_prefix> -p <pattern_file>
```

The command above shows to run queries using the document array profiles data-structure. It takes in the same output prefix
specified during the build process, along with FASTA file filled with reads that you want to query the input text for.

**Example of the output file:**
```sh
>pattern_0
[0,4] {0,1} 
>pattern_1
[0,2] {0,2} 
>pattern_2
[0,6] {0} 
>pattern_3
[0,4] {2} 
>pattern_4
[0,2] {0,1} 
>pattern_5
[2,5] {0} [1,1] {0,1,2} [0,0] {0,1} 
```
The `run` will create an output file ending in `.listings`, an example is shown above.
Each line in the output file represents a single read in the pattern file. The indexes in the square brackets
represent the substring in the read that were being search so `[0, 4]` represents a 5-bp string starting at position
 0 and including position 4. Then, the values in the braces represent the classes that the exact match
 occurs in, so `{0, 1}` means the 5-bp string occurs in document 0 and 1.
