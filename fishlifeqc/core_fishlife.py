#!/usr/bin/env python

# import sys
import argparse
from fishlifeqc.pairedblast import Pairedblast
from fishlifeqc.missingdata import Missingdata
from fishlifeqc.boldsearch  import Boldesearch
from fishlifeqc.t_like import TreeExplore
from fishlifeqc.bl import BLCorrelations
# from fishlifeqc.genetrees   import Raxml

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

    * horizontal and vertical trimming (default):

        $ fishlifeqc mdata [exons]

    * vertical trimming by codons: 

        $ fishlifeqc mdata [exons] --codon_aware

        note: stop codons are also counted along the sequence.
              if a particular sequence has more than one stop
              codon, the entire sequence is cut off from the 
              alignment 

    * vertical trimming by triplets and using mitochondrial 
      vertebrate lib of stop codons: 

        $ fishlifeqc mdata [exons] --codon_aware --mt

    * vertical trimming only:

        $ fishlifeqc mdata [exons] -v
""")

missingdata.add_argument('sequences', 
                    metavar='exons',
                    nargs="+",
                    type=str,
                    help='''File names with sequences''')
missingdata.add_argument('-e','--edges', 
                    metavar="",
                    type = float,
                    default = 0.6,
                    help='''[Optional] Set the vertical trimming threshold. 
                            Sequence columns at edges with more than 
                            this value are removed [Default: %s]''' % 0.6)
missingdata.add_argument('-c','--coverage', 
                    metavar="",
                    type = float,
                    default = 0.5,
                    help='''[Optional] Set the horizontal trimming threshold.
                             Sequences with more than this value are 
                             removed [Default: %s]''' % 0.5)
missingdata.add_argument('-t','--codon_aware',
                    action="store_true",
                    help='[Optional] If selected, trimming is done by codons')
missingdata.add_argument('-m','--mt',
                    action="store_true",
                    help='''[Optional] If `-t` is selected and this option are selected,
                     the mitochondrial vertebrate lib of stop codonds is used''')
missingdata.add_argument('-v','--verticaltrim',
                    action="store_false",
                    help='[Optional] If selected, only trim at adges is applied')
missingdata.add_argument('-q', '--quiet',
                    action='store_true',
                    help='[Optional] If selected, running messages are suppressed')
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
value is the query's group.  But, only if other group is detected
 (i.e. with mismatch sequence), this one is reported. 
Furthermore, a new file is created with sequences without
mismatches and any gaps produced by sequence elimination are 
also closed


Example:

    * Reciprocal blastn with a threshold of 99.9%:

        $ fishlifeqc rblast [exons] -t [taxonomy file] -i 99.9

            note 1: The taxnomy file is CSV-formated and must 
                    contain the following:

                    names               group              
                    [sequence header 1],[group of header 1]
                    [sequence header 2],[group of header 2]
                    [sequence header 3],[group of header 3]
                     ...                 ...            

            note 2: By default, this command will create a file 
                    with mismatches called `{0}`. The format of 
                    this file is CSV-formated and will 
                    contain the following:

                    exon ,sample             ,group
                    file1,[sequence header 1],[group of header 1] 
                    file1,[sequence header 2],[group of header 2]
                    file2,[sequence header 1],[group of header 1] 
                    file2,[sequence header 2],[group of header 2]
                    ...  ,      ...          ,     ...


    * Close gaps by codons:

        $ fishlifeqc rblast [exons] -t [taxonomy file] --codon_aware

""".format(PB_OUTPUTFILENAME))

pairedblast.add_argument('sequences', 
                    metavar='exons',
                    nargs="+",
                    type=str,
                    help='File names with sequences. If these are aligned, an unalignment process is performed')
pairedblast.add_argument('-t','--taxonomy',
                    metavar="",
                    default = None,
                    required= True,
                    help='Taxonomy file. Format in csv: "[sequence name],[group]"')
pairedblast.add_argument('-i','--identity', 
                    metavar="",
                    type = float,
                    default = PB_THRESOLD,
                    help='[Optional] Minimum identity values to perform each reciprocal blastn [Default: %s]' % PB_THRESOLD)
pairedblast.add_argument('-c','--codon_aware',
                    action="store_true",
                    help='[Optional] If selected, gaps produced by sequence elimination are closed by codons')
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
                                    help = "Match sequences against the BOLD database",
                                    formatter_class = argparse.RawDescriptionHelpFormatter, 
                                    description="""


                Wrapper of both BOLD and NCBI APIs for species identifications
                                    from DNA sequences
- Hosts:
    BOLD: http://www.boldsystems.org/index.php/Ids_xml
    NCBI: https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi

Example:

    $ fishlifeqc bold [sequence] -t [taxonomy file]

        The taxnomy file is CSV-formated and must contain the following:

            names               spps
            [sequence header 1],[species names of header 1]
            [sequence header 2],[species names of header 2]
            [sequence header 3],[species names of header 3]
             ...                 ...                    

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


# tlike ------------------------------------------------------

tlike = subparsers.add_parser('tlike',
                                    help = "Find T-like clades in trees",
                                    formatter_class = argparse.RawDescriptionHelpFormatter, 
                                    description="""


                T-like finder in trees

T-like topology is formed when two or more taxa have near-zero terminal
branch lengths in the same terminal clade:

                     |A
                -----|
                     |B

This near-zero value can be defined by the option `--min_len`.

Since many reconstruction methods do not consider polytomies, t-like clades
can be inside another t-like clade. Then, this command also collapse internal
branches in function of either length or support values to look for wider 
range of t-like clades.

By default, the gene trees are rooted by using the midpoint distance rooting 
method, but you can also use specific outgroups (see options).

Examples:

    * Standard usage:

        $ fishlifeqc tlike [tree files] -t [taxonomy file]

            Taxonomy file is CSV-formated and must contain the following
            struncture:

                # names       group 0                 group 1                 ... 
                [tip name 1],[group 0 of tip name 1],[group 1 of tip name 1], ...
                [tip name 2],[group 0 of tip name 2],[group 1 of tip name 2], ...
                [tip name 3],[group 0 of tip name 3],[group 1 of tip name 3], ...
                ...          ...                    ,...                  , ...

            Where group 0 can represent species names, group 1 genus, etc.
            NA values are empty spaces. 

    * Collapse internal branches by branch lengths:

        $ fishlifeqc tlike [tree files] -t [taxonomy file] -l 

    * Specify outgroups:

        $ fishlifeqc tlike [tree files] -t [taxonomy file] -g [outgroup file]

            Outgroup file is simple list with tip names in one column. 
            However, it can also contain a priorities in outgroups using two
            columns:

                 # priority, names
                 0,[tip name 1]
                 0,[tip name 2]
                 1,[tip name 3]
                 1,[tip name 4]
                 ...

            Where the first column represent group priority to outgroup gene trees.
            The lesser this number, higher the priority to choose those tip names
            for rooting gene trees. For example, if any tip of priority '0' is found
            at a gene tree, it is chosen tips of priority '1' to root the gene tree.
""")
tlike.add_argument('treefiles',
                    nargs="+",
                    help='Filenames')
tlike.add_argument('-t','--taxonomyfile',
                    metavar="",
                    type= str,
                    help='Taxonomy file [Default: None]') 
tlike.add_argument('-g','--outgroup',
                    metavar="",
                    type= str,
                    help='[Optional] Outgroup taxa [Default: None]')
tlike.add_argument('-f','--format',
                    metavar="",
                    type= str,
                    default= "newick",
                    help='[Optional] Tree format [Default: newick]') 
tlike.add_argument('-l','--coll_bylen',
                    action="store_true",
                    help='''[Optional] If selected, collapse internal branches by length''')
tlike.add_argument('-u','--ucoll_bysupp',
                    action="store_false",
                    help='''[Optional] If selected, internal branches are not collapsed by support value''')
tlike.add_argument('-L', '--min_len',
                    metavar = "",
                    type    = float,
                    default = 0.000001,
                    help    = '''[Optional] minimun branch length to collapse 
                    internal branch. This value is also used to 
                    define T-like clades [Default: 0.000001]''')
tlike.add_argument('-S', '--min_supp',
                    metavar = "",
                    type    = float,
                    default = 0,
                    help    = '[Optional] minimun support value to collapse internal branch [Default: 0]')
tlike.add_argument('-o','--outfile',
                    metavar="",
                    type= str,
                    default= "t_like.csv",
                    help='[Optional] Out filename [Default: t_like.csv]')
tlike.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default: 1]')   


# bl ------------------------------------------------------

bl = subparsers.add_parser('bl',
                            help = "Branch length ratios and correlations",
                            formatter_class = argparse.RawDescriptionHelpFormatter, 
                            description="""


                Branch length rations and pearson correlations

Both metrics are obtained by comparing a constrained
trees with a pruned reference tree

Example:

    * Standar usage:

        $ fishlifeqc bl [sequences] -t [reference tree]
""")

bl.add_argument('sequences',
                nargs="+",
                type=str,
                help='Filenames')
bl.add_argument('-t','--reference',
                metavar="n",
                type= str,
                default=None,
                required= True,
                help='Reference tree in newick format [Default: None]')
bl.add_argument('-r', '--max_ratio',
                 metavar = "",
                 type    = float,
                 default = 5.0,
                 help    = '[Optional] Constrained/pruned tree maximum ratio allowed [Default: 5]')
bl.add_argument('-p', '--min_pearson',
                 metavar = "",
                 type    = float,
                 default = 0.5,
                 help    = '[Optional] Minimum pearson correlation allowed [Default: 0.5]')
bl.add_argument('-e','--evomol',
                 metavar="",
                 type= str,
                 default = 'GTRGAMMA',
                 help='[Optional] RAxML evol. model for constrained tree inference [Default: None]')
bl.add_argument('-i', '--iterations',
                metavar = "",
                type    = int,
                default = 1,
                help    = '[Optional] Number of iterations for MLEs [Default: 1]')
bl.add_argument('-f','--prefix',
                metavar="",
                type= str,
                default= "BL_",
                help='[Optional] Prefix for output files [Default: BL_]')
bl.add_argument('-n', '--threads',
                metavar = "",
                type    = int,
                default = 1,
                help    = '[Optional] number of cpus [Default: 1]')

def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "mdata":
        # print(wholeargs)
        Missingdata(
                fastas       = wholeargs.sequences, 
                htrim        = wholeargs.coverage, 
                vtrim        = wholeargs.edges, 
                outputsuffix = wholeargs.out, 
                horizontal   = wholeargs.verticaltrim, # default true
                codon_aware  = wholeargs.codon_aware, # default false
                mtlib        = wholeargs.mt, # default true
                quiet        = wholeargs.quiet,
                threads      = wholeargs.threads
            ).run()

    elif wholeargs.subcommand == "rblast":

        Pairedblast(
            sequences = wholeargs.sequences,
            taxonomy  = wholeargs.taxonomy,
            threads   = wholeargs.threads,
            outname   = wholeargs.out,
            threshold = wholeargs.identity
        ).run()

    elif wholeargs.subcommand == "bold":
        
        Boldesearch(
            sequence           = wholeargs.sequence,
            bolddatabase       = wholeargs.bold_db,
            make_blast         = wholeargs.ncbi,
            quiet              = wholeargs.quiet,
            taxonomyfile       = wholeargs.taxonomy,
            removeintermediate = wholeargs.keep,
            threshold          = wholeargs.threshold,
            outfile            = wholeargs.out
        ).id_engine()

    elif wholeargs.subcommand == "tlike":
        # print(wholeargs)
        
        TreeExplore(
            treefiles      = wholeargs.treefiles,
            schema         = wholeargs.format,
            collpasebylen  = wholeargs.coll_bylen, # default: false
            minlen         = wholeargs.min_len,
            collpasebysupp = wholeargs.ucoll_bysupp, # default: true (counter pattern)
            minsupp        = wholeargs.min_supp,
            taxnomyfile    = wholeargs.taxonomyfile,
            outgroup       = wholeargs.outgroup,
            # suffix         = wholeargs.suffix, 
            threads        = wholeargs.threads,
            outfilename    = wholeargs.outfile
        ).find_Tlikes()

    elif wholeargs.subcommand == 'bl':

        BLCorrelations(
            species_tree_file = wholeargs.reference,
            sequences         = wholeargs.sequences,
            ratio_threshold   = wholeargs.max_ratio,
            pearson_threshold = wholeargs.min_pearson,
            iterations        = wholeargs.iterations,
            evomodel          = wholeargs.evomodel,
            prefix            = wholeargs.prefix,
            threads           = wholeargs.threads
        ).BrLengths()

    # elif wholeargs.subcommand == "raxmltree":
    #     # print(wholeargs) 
    #     Raxml(
    #         alignments     = wholeargs.exonfiles,
    #         concatenate    = wholeargs.concatenate, # false
    #         name_concate   = wholeargs.matrixname,
    #         evomodel       = wholeargs.model,
    #         bootstrap      = wholeargs.bootstrap,
    #         threads        = wholeargs.threads,
    #         raxml_failures = wholeargs.raxmlfailures,
    #     ).run()

if __name__ == "__main__":
    main()
