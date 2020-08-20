#!/usr/bin/env python3

import sys
import argparse
from fishlifeqc.pairedblast import Pairedblast
from fishlifeqc.missingdata import Missingdata
from fishlifeqc.boldsearch import Boldesearch


PB_MAKEBLASTFAILURE = "failed_to_makeblastdb.txt"
PB_OUTPUTFILENAME   = "mismatch_pairedblastn.txt"
PB_THRESOLD         = 95.0
BS_BOLD_DB          = 'COX1_SPECIES_PUBLIC'
BS_THRESHOLD        = 0.98

parser = argparse.ArgumentParser( formatter_class = argparse.RawDescriptionHelpFormatter, 
                                      description = '''
                                 Quality Control Steps
                                      ''')

subparsers = parser.add_subparsers(help='', dest='subcommand')


# mdata       ------------------------------------------------------
missingdata = subparsers.add_parser('mdata',
                                    help = "Trim sequences in function of gap ocurrences",
                                    formatter_class = argparse.RawDescriptionHelpFormatter, 
                                    description="""
                            Missing data

Examples:

    * horizontal trimming with a threshold of 0.5 (default):

        $ fishlifeqc mdata [exons] -c 0.5

    * vertical trimming with a threshold value of 0.6 (default):

        $ fishlifeqc mdata [exons] -v -e 0.6

    * vertical trimming by triplets: 

        $ fishlifeqc mdata [exons] -v --codon_aware

        note: stop codons are also counted along the sequence.
              if a particular sequence has more than one stop
              codon, the entire sequence is cut off from the 
              alignment 

    * vertical trimming by triplets and using mitochondrial 
      vertebrate lib of stop codons: 

        $ fishlifeqc mdata [exons] -v --codon_aware --mt
""")

missingdata.add_argument('sequences', 
                    metavar='',
                    nargs="+",
                    type=str,
                    help='''File names with sequences''')
missingdata.add_argument('-c','--coverage', 
                    metavar="",
                    type = float,
                    default = 0.5,
                    help='''[Optional] Set the horizontal trimming threshold.
                             Sequences with more than this value are 
                             removed [Default: %s]''' % 0.5)
missingdata.add_argument('-v','--verticaltrim',
                    action="store_true",
                    help='[Optional] if selected, trim at adges is applied')
missingdata.add_argument('-e','--edges', 
                    metavar="",
                    type = float,
                    default = 0.6,
                    help='''[Optional] If `-v` is selected, set the vertical trimming threshold. 
                            Sequence columns at edges with more than 
                            this value are removed [Default: %s]''' % 0.6)
missingdata.add_argument('-t','--codon_aware',
                    action="store_true",
                    help='[Optional] If `-v` and this option are selected, trimming is done by triplets')
missingdata.add_argument('-m','--mt',
                    action="store_true",
                    help='''[Optional] If `-a` is selected and this option are selected,
                     the mitochondrial vertebrate lib of stop codonds is used''')
missingdata.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')                        
missingdata.add_argument('-o','--out', 
                    metavar="",
                    default = "_trimmed",
                    type= str,
                    help='[Optional] Suffix for outfile name [Default: %s]' % "_trimmed")

# mdata       ------------------------------------------------------

# pairedblast ------------------------------------------------------
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

                                        """ % PB_OUTPUTFILENAME)

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
                    default = PB_THRESOLD,
                    help='[Optional] Minimum identity values to perform each reciprocal blastn [Default: %s]' % PB_THRESOLD)
pairedblast.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')                        
pairedblast.add_argument('-o','--out', 
                    metavar="",
                    default = PB_OUTPUTFILENAME,
                    type= str,
                    help='[Optional] output file [Default: %s]' % PB_OUTPUTFILENAME)
# pairedblast ------------------------------------------------------

# bold search ------------------------------------------------------
boldsearch = subparsers.add_parser('bold',
                                    help = "Look for matches between sequences and the BOLD database",
                                    formatter_class = argparse.RawDescriptionHelpFormatter, 
                                    description="""


                Wrapper of both BOLD and NCBI APIs for species identifications
                                    from DNA sequences
- Hosts:
    BOLD: http://www.boldsystems.org/index.php/Ids_xml
    NCBI: https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi

- Example:
    fishlifeqc bold [sequence] -t [taxonomy file]
""")
boldsearch.add_argument('sequence', 
                        metavar='',
                        type=str,
                        help='''File name with the COI sequences. 
                                If these are aligned, an 
                                unalignment process is performed''')
boldsearch.add_argument('-t','--taxonomy',
                        metavar="",
                        default = None,
                        required= True,
                        help='Taxonomy file. Format in csv: [sequence name],[group],[species name]')
boldsearch.add_argument('-v','--threshold',
                        type = float,
                        metavar="",
                        action='store',
                        default=BS_THRESHOLD,
                        help='Minimum similarity allowed for best matched species [Default = %s]' % BS_THRESHOLD)
boldsearch.add_argument('-b','--bold_db', 
                        metavar="",
                        type  = str,
                        default = BS_BOLD_DB,
                        action='store',
                        help='''BOLD database. There are four available: 
                                COX1,
                                COX1_SPECIES,
                                COX1_SPECIES_PUBLIC,
                                COX1_L640bp
                                [Default = %s]''' % BS_BOLD_DB )
boldsearch.add_argument('-n','--ncbi',
                        action='store_true',
                        help=' If selected, BLASTn is used to identify species')
boldsearch.add_argument('-k','--keep',
                        action='store_false',
                        help=' If selected, intermideate files are kept')
boldsearch.add_argument('-q', '--quiet',
                        action='store_true',
                        help=' If selected, suppress running messages')
boldsearch.add_argument('-o','--out', 
                        metavar="",
                        type = str,
                        default=None,
                        action='store',
                        help='Output name [Default = `sequence` + _bold ]' )
# bold search ------------------------------------------------------

def main():

    wholeargs = parser.parse_args()
    if wholeargs.subcommand == "mdata":
        # print(wholeargs)
        Missingdata(
                fastas = wholeargs.sequences, 
                htrim = wholeargs.coverage, 
                vtrim = wholeargs.edges, 
                outputsuffix = wholeargs.out, 
                trimedges = wholeargs.verticaltrim, # default false
                codon_aware = wholeargs.codon_aware, # default false
                mtlib = wholeargs.mt, # default true
                threads = wholeargs.threads
            ).run()

    elif wholeargs.subcommand == "rblast":

        pblast = Pairedblast(
                        sequences = wholeargs.sequences,
                        taxonomy  = wholeargs.taxonomy,
                        threads   = wholeargs.threads,
                        threshold = wholeargs.identity
                        )

        pblast.checktaxonomyfile()

        passed, failed = pblast.iteratedbproduction()

        if failed:
            with open(PB_MAKEBLASTFAILURE, "w") as f:
                for i in failed:
                    f.write(i + "\n")

        if passed:

            outliers = pblast.iterateblastn(passed)

            if outliers:
                
                with open(wholeargs.out, 'w') as f:

                    f.write("exon,sample,group\n")

                    for exon,spps,group in outliers:
                        f.write(  "%s,%s,%s\n" %  (exon, spps, group) )

    elif wholeargs.subcommand == "bold":
        
        Boldesearch(
            sequence = wholeargs.sequence,
            bolddatabase = wholeargs.bold_db,
            make_blast = wholeargs.ncbi,
            quiet = wholeargs.quiet,
            taxonomyfile = wholeargs.taxonomy,
            removeintermediate = wholeargs.keep,
            threshold = wholeargs.threshold,
            outfile = wholeargs.out
        ).id_engine()

if __name__ == "__main__":
    main()