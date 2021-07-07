# fishlifeqc

Features:

- [x] quality control for exon alignments

Software requierements:

* pip

## Installation

Main version by using `pip`:
```Bash
pip install fishlifeqc
```

Development version by using `git` and `pip` (optional):
```Bash
git clone https://github.com/Ulises-Rosas/fishlifeqc.git 
cd fishlifeqc
python3 -m pip install .
```

## Usage

Main command:
```Bash
fishlifeqc -h
```

```
usage: fishlifeqc [-h] {mdata,rblast,bold,tlike,bl,para,srh} ...

                                 Quality Control Steps
                                      

positional arguments:
  {mdata,rblast,bold,tlike,bl,para,srh}
    mdata               Trim sequences in function of gap ocurrences
    rblast              Reciprocal blastn comparing taxonomical groups
    bold                Match sequences against the BOLD database
    tlike               Find T-like clades in trees
    bl                  Branch length ratios and correlations
    para                Test paraphyly with AU tests
    srh                 Test stationarity, reversibility and homogeneity

optional arguments:
  -h, --help            show this help message and exit
```

Utilities command:
```Bash
qcutil -h
```

```
usage: qcutil [-h] {qstats,fstats,itt,rf,delete,knockdown} ...

                                 Utilities from fishlifeqc
                                      

positional arguments:
  {qstats,fstats,itt,rf,delete,knockdown}
    qstats              quick summary of alignment information
    fstats              full summary of both alignment and tree information
    itt                 Incongruence through time
    rf                  Robinson-Foulds distances
    delete              Delete specific sequences from files
    knockdown           Replace specific codon positions with 'N's

optional arguments:
  -h, --help            show this help message and exit
```

A description of summary statistics currently available for `qcutil` can be found [here](https://github.com/Ulises-Rosas/fishlifeqc/blob/master/qcutil/var_names.md).


Yet unintegrated scripts into `qcutil`:

* `merge.py` : join two set of exons by using a map file
* `concatenate.py` : concatenate exons
* `codon_partition.py` : generate a condon partition file
* `cons_trees.py` : for each exon, constrain a reference tree by exon taxa
* `collapse.py` : collapse edges in function of support values
