#!/usr/bin/env python3

import os
import re
import argparse
from fishlifeqc.bl import BLCorrelations

def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

            For each exon file, a constrained tree is generated
            by intersecting tree taxa and the exon file taxa

    * Usage:

        $ cons_trees.py [exon files] -t [tree]

""")
    parser.add_argument('filenames',
                        nargs="+",
                        help='exons names')
    parser.add_argument('-t','--tree_file',
                        metavar="",                                       
                        required=True,
                        type= str,
                        default= None,
                        help='Tree file [Default: None]')   
    parser.add_argument('-b','--bl',
                        action="store_true",
                        help='''[Optional] If selected, only constrained trees with branch lengths 
                                are generated''')
    parser.add_argument('-s', '--suffix',
                        metavar = "",
                        type    = str,
                        default = '_constr.tree',
                        help    = '[Optional] Suffix for output files [Default: _constr.tree]')
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default: 1]')

    args = parser.parse_args()
    return args

def main():
    args = getOpts()
    # print(args)
    BLCorrelations(
        species_tree_file = args.tree_file,
        sequences = args.filenames,
        with_bl = args.bl,
        suffix  = args.suffix,
        threads = args.threads
    ).get_constraints()
    
if __name__ == "__main__":
    main()
