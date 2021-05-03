#!/usr/bin/env python

import argparse
from qcutil.summarize import Stats,Incongruence

ITT_OUTFILENAME = "support_tt.tsv"
RF_OUTFILENAME = "RF_distance.tsv"

parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                                 Utilities from fishlifeqc
                                      ''')
subparsers = parser.add_subparsers(help='', dest='subcommand')


# stats --------------------------------------------------------------------------
describe = subparsers.add_parser('stats', help = "Summarize alignment information",
                              formatter_class = argparse.RawDescriptionHelpFormatter, 
                              description="""
                              
                    Summarize alignment information

Examples:

    * Standard usage:

        $ %(prog)s [alignment files]
""")

describe.add_argument('filenames',
                      metavar = 'file',
                      nargs="+",
                      help='Filenames')
describe.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')                        
# stats --------------------------------------------------------------------------

# itt --------------------------------------------------------------------------
itt = subparsers.add_parser('itt', help = "Incongruence through time",
                              formatter_class = argparse.RawDescriptionHelpFormatter, 
                              description="""
                              
        Examination of support and conflict of gene trees accross nodes

Examples:

    * Standard usage:

        $ %(prog)s [gene trees] -r [reference tree]
""")
itt.add_argument('filenames',
                  metavar = 'file',
                  nargs="+",
                  help='Filenames')
itt.add_argument('-r', '--ref_tree',
                  metavar = '',
                  required=True,
                  help='A dated tree in newick format')
itt.add_argument('-o', '--out',
                  metavar = '',
                  type=str,
                  default=ITT_OUTFILENAME,
                  help='[Optional] output file [Default: %s]' % ITT_OUTFILENAME)
itt.add_argument('-n', '--threads',
                  metavar = "",
                  type    = int,
                  default = 1,
                  help    = '[Optional] number of cpus [Default = 1]')
# itt --------------------------------------------------------------------------


# RF --------------------------------------------------------------------------
rf = subparsers.add_parser('rf', help = "Robinson-Foulds distances",
                              formatter_class = argparse.RawDescriptionHelpFormatter, 
                              description="""
                              
             Robinson-Foulds distances

Examples:

    * Standard usage:

        $ %(prog)s [gene trees] -r [reference tree]
""")
rf.add_argument('filenames',
                 metavar = 'file',
                 nargs="+",
                 help='Filenames')
rf.add_argument('-r', '--ref_tree',
                 metavar = '',
                 required=True,
                 type=str,
                 help='A dated tree in newick format')
rf.add_argument('-w','--weighted',
                 action='store_true',
                 help=' If selected, weighted Robison-Foulds distances are calculated')
rf.add_argument('-o', '--out',
                 type=str,
                 metavar = '',
                 default=RF_OUTFILENAME,
                 help='[Optional] output file [Default: %s]' % RF_OUTFILENAME)
rf.add_argument('-n', '--threads',
                 metavar = "",
                 type    = int,
                 default = 1,
                 help    = '[Optional] number of cpus [Default = 1]')
# RF --------------------------------------------------------------------------

def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "stats":

        Stats(
            fastas= wholeargs.filenames,
            threads= wholeargs.threads
        ).run()

    elif wholeargs.subcommand == "itt":

        Incongruence(
            gene_trees     = wholeargs.filenames,
            reference_tree = wholeargs.ref_tree,
            quiet          = False,
            threads        = wholeargs.threads

        ).support_tt( outfile = wholeargs.out )

    elif wholeargs.subcommand == "rf":
        # print(wholeargs)
        Incongruence(
            gene_trees     = wholeargs.filenames,
            reference_tree = wholeargs.ref_tree,
            weighted_rf    = wholeargs.weighted,
            quiet          = False,
            threads        = wholeargs.threads

        ).get_rf( outfile = wholeargs.out )


if __name__ == "__main__":
    main()
