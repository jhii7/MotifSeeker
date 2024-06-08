# MotifSeeker

This was created as a final project for Professor Melissa Gymrek's CSE 185 class at UCSD. MotifSeeker is a python command-line tool which aims to identify enriched motifs in the genome when provided reference genome and peak files. It aims to emulate some limited functionality of the HOMER python package.

## Installation

`MotifSeeker` requires the following python libraries to be installed:
- numpy
- argparse
- bed_reader
- scipy
- seqlogo
- biopython

They can be installed using `pip`:

```pip install argparse numpy bed-reader scipy seqlogo biopython```

### Install instructions

To install our tool, you can run the following commands:

```
git clone https://github.com/jhii7/MotifSeeker.git
cd MotifSeeker
python setup.py install
```
Make sure you are in the `~/MotifSeeker` directory when you run the last command.

If the install was successful, typing `motifseeker --help` should show a useful message.

### Usage instructions

`MotifSeeker` uses the following input files from the command-line arguments:
- `.bed file` containing the peak data
- `.fa file` containing the reference genome

### Basic usage

The basic usage of `motifseeker` is:
```
motifseeker <.bed file> <.fa ref genome file> [other options]
```

To run `motifseeker` on a small test example (using files in this repo):
```
motifseeker example_files/smallbed.bed example_files/test.fa
```
Make sure you are in the `~/MotifSeeker` directory when you run the above command.

### Additional command line option

`-p PVAL`, `--pval PVAL`: p-value threshold for significant enrichment. Default: 0.05

### Input file format requirements

The BED file should have a minimum of 3 tab separated columns (additional columns will be ignored)
(Read the [HOMER documentation](http://homer.ucsd.edu/homer/ngs/peakMotifs.html#:~:text=The%20findMotifsGenome.pl%20program%20is,the%20enrichment%20of%20known%20motifs.) for more on acceptable input files)
- Column1: chromosome number in the format `chr*` where `*` is the chromosome number
- Column2: starting position
- Column3: ending position

The fasta file should have the following format
```
>chr[name]
[chromosome sequence]
```

To refer to the standard motifs we use, refer to the file under `motifs` > `custom.motifs`

## Contributors

This repository was created by Tusha Karnani, Justin Hii, and Numa Yazadi. We would like to thank Dr. Gymrek for her [CSE 185 demo project](https://github.com/gymreklab/cse185-demo-project).
