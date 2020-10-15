# fishlifeqc

Features:

- [x] quality control for exon alignments

Software requierements:

* Python 3
* Git

## Installation

Get the development version:

```Bash
git clone https://github.com/Ulises-Rosas/fishlifeqc.git
cd fishlifeqc
python3 setup.py install
python3 setup.py install_data
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
