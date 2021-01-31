# fishlifeqc

Features:

- [x] quality control for exon alignments

Software requierements:

* pip
* git

## Installation

Get the development version:

```Bash
git clone https://github.com/Ulises-Rosas/fishlifeqc.git && cd fishlifeqc
pip install .
```

## Usage

```Bash
fishlifeqc -h
```

```
usage: fishlifeqc [-h] {mdata,rblast,bold,tlike,bl} ...

                                 Quality Control Steps
                                      

positional arguments:
  {mdata,rblast,bold,tlike,bl}
    mdata               Trim sequences in function of gap ocurrences
    rblast              Reciprocal blastn comparing taxonomical groups
    bold                Match sequences against the BOLD database
    tlike               Find T-like clades in trees
    bl                  Branch length ratios and correlations

optional arguments:
  -h, --help            show this help message and exit

```
