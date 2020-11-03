# BSc project (computer science - 2020/2021)

Software Design Project is a course held at University of Zagreb, Faculty of Electrical Engineering and Computing in the fifth semester of the undergraduate study. The main focus is to promote cooperation between students while solving specific problems. Under the supervision of professor Mile Šikić and assistant professor Krešimir Križanović, students will get familiar with C++, basics of compilation methods, version control, unit tests and continuous integration, and will be introduced to algorithms used in bioinformatics. This project will be executed in several steps, each with defined outcomes which are required for succeeding steps. Instructions and guidelines for the project will be written in this README file which will be updated through the semester.

## Preliminaries

Students are required to get through the following tutorials: [C++](http://www.cplusplus.com/doc/tutorial/), [GitHub](http://rogerdudler.github.io/git-guide/), [CMake](https://cmake.org/cmake/help/latest/guide/tutorial/index.html), [Googletest](https://github.com/google/googletest/blob/master/googletest/docs/primer.md) and [TravisCI](https://docs.travis-ci.com/user/getting-started/). While writing C++ code, we advise to follow the [Google C++ style guide](https://google.github.io/styleguide/cppguide.html).

Students will be assigned to one of five teams which are **blonde**, **blue**, **brown**, **orange**, **pink** and **white**. Each team will have a separate branch and only team members will have write permission. Although, additional branches can be created if needed, but should have names starting with the team name (e.g. `blonde_feature_one`).

## Objective

At the end of the project, students will have implemented several libraries that enable placement and description of similarity between a large amount of substrings (of various sizes) and a much longer string from which they all originate. The goal is to create a single program out of those libraries in order to map long erroneous fragments obtained via third generation sequencing technologies to a reference genome, which has various use cases in bioinformatics. A visual example can be seen bellow.

![](misc/sample_mappings.png)

## Setup

Each team's main branch should be up to date with this README. The project setup consists of creating a C++ program and naming it in form of `<team name>_mapper` (e.g. `blonde_mapper`) with `cmake`. It has to accept two files as floating arguments and enable options `-h` (`--help`) and `--version`, which are used for the help and version messages (following [SemVer](https://semver.org/)), respectively. Suggested argument parser to include is `optarg`, but this feature can also be implemented independently.

The first file will contain a reference genome in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format, while the second file will contain a set of fragments in either FASTA or [FASTQ](https://en.wikipedia.org/wiki/FASTQ_format) format. The files need to be parsed and stored in memory, and some statistics have to be outputted to `stderr`, which includes names of sequences in the reference file and their lengths, number of sequences in the fragments file, their average length, N50 length, minimal and maximal length, etc. There is no need to implement a parser, you can add [bioparser](https://github.com/rvaser/bioparser) to the project as a submodule via `git` and integrate it with `cmake`. It supports several bioinformatics formats where files can also be compressed with `gzip`.

Sample program runs, after the setup step is completed, can be seen bellow:

```bash
blonde_mapper GCF_000005845.2_ASM584v2_genomic.fna MAP006-1_2D_pass.fasta
<basic statistics of input files>
```

```bash
blonde_mapper -h
<appropriate message describing supported arguments>
```

```bash
blonde_mapper --version
v0.1.0
```

## Data

The first version of the mapper will be tested on an Oxford Nanopore Technologies data set obtained by sequencing the Escherichia coli K-12 substr. MG1655 genome. The data set is available from Loman Labs [here](http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/) (download one or both of MAP-006-1 and MAP-006-2 FASTA files), while the reference genome is available from NCBI [here](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz).
