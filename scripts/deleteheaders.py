
#!/usr/bin/env python3

import argparse
from fishlifeqc.delete import Deletion

def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

            Delete specific sequences from multiple files

    * Usage:

        $ deleteheaders.py [exon files] -c [control file]

            Where the control file is a simple list of 
            sequences to be deleted in plain text

    * Specify number of threads:

        $ deleteheaders.py [exon files] -c [control file] -n 5


""")
    parser.add_argument('filenames',
                        nargs="+",
                        help='Filenames')
    parser.add_argument('-c','--controlfile',
                        metavar="",
                        default = None,
                        required= True,
                        help='Control file with the list of species')
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')                        
    args = parser.parse_args()
    return args

def main():
    args = getOpts()

    Deletion(
        sequences   = args.filenames,
        controlfile =  args.controlfile,
        filetype    =  'list',
        threads     = args.threads 
    ).headers()

if __name__ == "__main__":
    main()