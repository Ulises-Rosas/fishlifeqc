#!/usr/bin/env python3

import sys
import argparse
from fishlifeqc.pairedblast import Pairedblast


parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                                 Quality Control Steps
                                      ''')

subparsers = parser.add_subparsers(help='', dest='subcommand')

# pairedblast ------------------------------------------------------

MAKEBLASTFAILURE = "failed_to_makeblastdb.txt"
OUTPUTFILENAME   = "mismatch_pairedblastn.txt"
THRESOLD         = 95.0

pairedblast = subparsers.add_parser('rblast',
                                    help = "Reciprocal blastn comparing taxonomical groups",
                                    formatter_class = argparse.RawDescriptionHelpFormatter, 
                                    description="""

                    Reciprocal blastn comparing taxonomical groups

                The expected group for each blastn with a given threshold 
                value is the query's group.  But, only if other group is detected,
                this one is reported at '%s' by default. 
                Filename can be changed with `-o` option. 
                See below for further details. 

Example:

    fishlifeqc rblast [exons] -t [taxonomy file]

                                        """ % OUTPUTFILENAME)

pairedblast.add_argument('sequences', 
                    metavar='',
                    nargs="+",
                    type=str,
                    help='File names with sequences. If these are aligned, an unalignment process is performed')
pairedblast.add_argument('-t','--taxonomy',
                    metavar="",
                    default = None,
                    required= True,
                    help='Taxonomy file. Format in csv: [sequence name],[group],[species name]')
pairedblast.add_argument('-i','--identity', 
                    metavar="",
                    type = float,
                    default = THRESOLD,
                    help='[Optional] Minimum identity values to perform each reciprocal blastn [Default: %s]' % THRESOLD)
pairedblast.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')                        
pairedblast.add_argument('-o','--out', 
                    metavar="",
                    default = OUTPUTFILENAME,
                    type= str,
                    help='[Optional] output file [Default: %s]' % OUTPUTFILENAME)
# pairedblast ------------------------------------------------------

def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "rblast":

        pblast = Pairedblast(
                        sequences = wholeargs.sequences,
                        taxonomy  = wholeargs.taxonomy,
                        threads   = wholeargs.threads,
                        threshold = wholeargs.identity
                        )

        pblast.checktaxonomyfile()

        passed, failed = pblast.iteratedbproduction()

        if failed:
            with open(MAKEBLASTFAILURE, "w") as f:
                for i in failed:
                    f.write(i + "\n")

        if passed:

            outliers = pblast.iterateblastn(passed)

            if outliers:
                
                with open(wholeargs.out, 'w') as f:

                    f.write("exon,sample,group\n")

                    for exon,spps,group in outliers:
                        f.write(  "%s,%s,%s\n" %  (exon, spps, group) )
        
if __name__ == "__main__":
    main()