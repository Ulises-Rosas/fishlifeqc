#!/usr/bin/env python3

import os
import argparse
from fishlifeqc.codons import Codonm


def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

            Get horizontal (i.e. at sequence level) and 
        vertical (i.e. at alignment colum level) base composition
                for different exon files

	* Example:

	    $ codomco.py [exon files]
""",
                                     epilog="")
    parser.add_argument('filenames',
                        nargs="+",
                        help='Filenames')
    parser.add_argument('-s','--suffix',
                        metavar="",
                        type= str,
                        default= "",
                        help='[Optional] Suffix for codom composition files [Default = ""]') 
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')                        
    args = parser.parse_args()
    return args

def main():
    args = getOpts()

    Codonm(
        threads = args.threads,
        alns    = args.filenames,
        suffix  = args.suffix
    ).run()

if __name__ == "__main__":
    main()
