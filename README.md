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
usage: fishlifeqc [-h] {mdata,rblast,bold} ...

                                 Quality Control Steps
                                      

positional arguments:
  {mdata,rblast,bold}
    mdata              Trim sequences in function of gap ocurrences
    rblast             Reciprocal blastn comparing taxonomical groups
    bold               Match sequences against the BOLD database

optional arguments:
  -h, --help           show this help message and exit

```
