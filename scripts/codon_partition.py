#!/usr/bin/env python3

import argparse
from multiprocessing import Pool
from fishlifeqc.utils import codon_partitions_nexus

def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

            Create codon partition in Nexus format

	* Example:

	    $ codon_partition.py [exon files]
            
            note: partition file name will have the same
                  as the entered alignment + ".nex" extension
""",
                                     epilog="")
    parser.add_argument('filenames',
                        nargs="+",
                        help='Filenames')
    # parser.add_argument('-s','--suffix',
    #                     metavar="",
    #                     type= str,
    #                     default= "",
    #                     help='[Optional] Suffix for codom composition files [Default = ""]') 
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')                        
    args = parser.parse_args()
    return args

def main():
    args = getOpts()

    with Pool(processes = args.threads) as p:
        matched  = [*p.map(codon_partitions_nexus, args.filenames)]

    no_partitions = list(filter(None, matched))

    if no_partitions:    
        with open('no_partitions.txt', 'w') as f:
            f.writelines(no_partitions)

if __name__ == "__main__":
    main()