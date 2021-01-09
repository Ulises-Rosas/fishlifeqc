
#!/usr/bin/env python3

import argparse
from fishlifeqc.delete import Deletion

def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

            Delete specific sequences from multiple files

    * Standard usage:

        $ deleteheaders.py -i [exon files] -c [control file]

            Where the control file is a simple list of 
            sequences to be deleted in plain text

    * Specify number of threads:

        $ deleteheaders.py -i [exon files] -c [control file] -n 5

    * Delete specific headers per exon:

        $ deleteheaders.py -e -c [control file] -n 5

            Where `control file` has the following format:
                [exon 1],[header 1]
                [exon 1],[header 2]
                [exon 2],[header 1]
                [exon 3],[header 1]
                ...     ,...
""")
    parser.add_argument('-i','--filenames',
                        metavar="",
                        nargs="+",
                        help='Filenames')
    parser.add_argument('-e','--exon_header',
                        action="store_true",
                        help='''[Optional] If selected, delete by control file only''')
    parser.add_argument('-c','--controlfile',
                        metavar="",
                        default = None,
                        required= True,
                        help='[Optional] Control file with the list of species')
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')
    parser.add_argument('-s','--suffix',
                        metavar="",
                        type= str,
                        default= "_listd",
                        help='[Optional] Suffix for outputs [Default: _listd]')
    args = parser.parse_args()
    return args

def main():
    args = getOpts()
    # print(args)
    if args.exon_header:
        Deletion(
            sequences   = None,
            controlfile = args.controlfile,
            filetype    = args.suffix,
            threads     = args.threads 
        ).header_exon()

    else:
        Deletion(
            sequences   = args.filenames,
            controlfile = args.controlfile,
            filetype    = args.suffix,
            threads     = args.threads 
        ).headers()

if __name__ == "__main__":
    main()