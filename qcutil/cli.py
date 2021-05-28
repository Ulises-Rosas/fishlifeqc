#!/usr/bin/env python

import argparse
from qcutil.summarize import Stats,Incongruence
from fishlifeqc.delete import Deletion


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
describe.add_argument('-a','--per_aln',
                        action='store_true',
                        help=' If selected, summary information is done per alignment')
describe.add_argument('-s','--per_seq',
                        action='store_true',
                        help=' If selected, summary information is done per sequence')
describe.add_argument('-p','--prefix', 
                        metavar="",
                        type = str,
                        default='stats',
                        help='prefix name for outputs [Default = stats ]' )
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

# delete --------------------------------------------------------------------------
delete = subparsers.add_parser('delete', 
                                help = "Delete specific sequences from files",
                                formatter_class = argparse.RawDescriptionHelpFormatter, 
                                description="""

            Delete specific sequences from multiple files

    * Standard usage:

        $ deleteheaders.py [exon files] -c [control file]

            Where the control file is a simple list of 
            sequences to be deleted in plain text

    * Specify number of threads:

        $ deleteheaders.py [exon files] -c [control file] -n 5

    * Delete specific headers per exon:

        $ deleteheaders.py -e -c [control file] -n 5

            Where `control file` has the following format:
                [exon 1],[header 1]
                [exon 1],[header 2]
                [exon 2],[header 1]
                [exon 3],[header 1]
                ...     ,...
""")
delete.add_argument('filenames',
                    metavar="",
                    type=str,
                    nargs="*",
                    help='Filenames')
delete.add_argument('-e','--exon_header',
                    action="store_true",
                    help='''[Optional] If selected, delete by control file only''')
delete.add_argument('-c','--controlfile',
                    metavar="",
                    default = None,
                    required= True,
                    help='[Optional] Control file with the list of species')
delete.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')
delete.add_argument('-s','--suffix',
                    metavar="",
                    type= str,
                    default= "_listd",
                    help='[Optional] Suffix for outputs [Default: _listd]')
# delete --------------------------------------------------------------------------

def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "stats":

        Stats(
            fastas      = wholeargs.filenames,
            align_based = wholeargs.per_aln,
            seq_based   = wholeargs.per_seq,
            prefix      = wholeargs.prefix,
            threads     = wholeargs.threads
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

    elif wholeargs.subcommand == "delete":

        if wholeargs.exon_header:
            Deletion(
                sequences   = None,
                controlfile = wholeargs.controlfile,
                filetype    = wholeargs.suffix,
                threads     = wholeargs.threads 
            ).header_exon()

        else:
            Deletion(
                sequences   = wholeargs.filenames,
                controlfile = wholeargs.controlfile,
                filetype    = wholeargs.suffix,
                threads     = wholeargs.threads 
            ).headers()

if __name__ == "__main__":
    main()
