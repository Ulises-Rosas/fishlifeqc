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

To look at the help page of `pairedblast`, you can use:

```Bash
fishlifeqc -h
```

```usage: fishlifeqc [-h] {mdata,rblast,bold,delete,raxmltree} ...

                                 Quality Control Steps
                                      

positional arguments:
  {mdata,rblast,bold,delete,raxmltree}
    mdata               Trim sequences in function of gap ocurrences
    rblast              Reciprocal blastn comparing taxonomical groups
    bold                Look for matches between sequences and the BOLD
                        database
    delete              Delete specific sequences from fasta files
    raxmltree           Get raxml trees for each exon

optional arguments:
  -h, --help            show this help message and exit
```
