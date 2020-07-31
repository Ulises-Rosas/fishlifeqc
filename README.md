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
pairedblast -h
```

```
usage: pairedblast [-h] -t  [-i] [-n] [-o]  [...]

                        Paired blastn comparing taxonomical groups

The expected group for each blastn with a given threshold value is the query's group. 
But, only if other group is detected, this one is reported at 'mismatch_pairedblastn.txt' by default. Filename
can be changed with `-o` option. See below for further details. 

positional arguments:
                    File names with sequences. If these are aligned, an
                    unalignment process is performed

optional arguments:
  -h, --help        show this help message and exit
  -t , --taxonomy   Taxonomy file. Format in csv: [sequence
                    name],[group],[species name]
  -i , --identity   [Optional] Minimum identity values to perform each
                    reciprocal blastn [Default: 95.0]
  -n , --threads    [Optional] number of cpus [Default = 1]
  -o , --out        [Optional] output file [Default:
                    mismatch_pairedblastn.txt]
```
