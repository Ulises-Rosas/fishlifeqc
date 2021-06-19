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
python3 -m pip install .
```

## Usage

Main command:
```Bash
fishlifeqc -h
```

```
usage: fishlifeqc [-h] {mdata,rblast,bold,tlike,bl,para} ...

                                 Quality Control Steps
                                      

positional arguments:
  {mdata,rblast,bold,tlike,bl,para}
    mdata               Trim sequences in function of gap ocurrences
    rblast              Reciprocal blastn comparing taxonomical groups
    bold                Match sequences against the BOLD database
    tlike               Find T-like clades in trees
    bl                  Branch length ratios and correlations
    para                Test paraphyly with AU tests

optional arguments:
  -h, --help            show this help message and exit

```

Utilities command:
```Bash
qcutil -h
```

```
usage: qcutil [-h] {qstats,fstats,itt,rf,delete} ...

                                 Utilities from fishlifeqc
                                      

positional arguments:
  {qstats,fstats,itt,rf,delete}
    qstats              quick summarize of alignment information
    fstats              full summarize of both alignment and tree information
    itt                 Incongruence through time
    rf                  Robinson-Foulds distances
    delete              Delete specific sequences from files

optional arguments:
  -h, --help            show this help message and exit
```
