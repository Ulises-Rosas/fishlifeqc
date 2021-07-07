#!/usr/bin/env python

# import sys
import argparse
from fishlifeqc.missingdata import Missingdata, STOP_CODON_TABLE

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

    * Codon aware trimming: 

        $ fishlifeqc mdata [exons] --codon_aware

        note: stop codons are also counted along the sequence.
              if a particular sequence has more than one stop
              codon, the entire sequence is cut off from the 
              alignment 

    * Trimming considering mitochondrial vertebrate lib of stop codons: 

        $ fishlifeqc mdata [exons] --codon_aware --stop_lib 2

        note: for a complete list of stop codon libraries run:
              $ fishlifeqc mdata -p
""")

missingdata.add_argument('sequences', 
                    metavar='exons',
                    nargs="*",
                    type=str,
                    help='''File names with sequences''')
missingdata.add_argument('-e','--edges', 
                    metavar="",
                    type = float,
                    default = 0.6,
                    help='''[Optional] Maximum allowed gap proportion per column at edges. 
                            Sequence columns at edges with more than 
                            this value are removed [Default: %s]''' % 0.6)
missingdata.add_argument('-i','--internal', 
                    metavar="",
                    type = float,
                    default = 0.6,
                    help='''[Optional] Maximum allowed gap proportion per internal column. 
                            Sequence internal columns with more than 
                            this value are removed [Default: %s]''' % 0.6)
missingdata.add_argument('-H','--coverage',
                    metavar="",
                    type = float,
                    default = 0.5,
                    help='''[Optional] Maximum allowed gap proportion per sequence.
                             Sequences with more than this value are 
                             removed [Default: %s]''' % 0.5)
missingdata.add_argument('-m','--min_seqs_per_aln',
                    metavar="",
                    type = int,
                    default = 4,
                    help='''[Optional] Minimum number of sequences per alignments. 
                            'qcutil stats' with the '-a' option might help to find this parameter [Default: %s]''' % 4)
missingdata.add_argument('-w','--min_alns_per_seq',
                    metavar="",
                    type = int,
                    default = 1,
                    help='''[Optional] Minimum number of alignments per sequence.
                            'qcutil stats' with the '-s' option might help to find this parameter [Default: %s]''' % 1)
missingdata.add_argument('-c','--codon_aware',
                    action="store_true",
                    help='[Optional] If selected, trimming is done by codons')
                    
missingdata.add_argument('-u','--unadjusted',
                    action="store_true",
                    help='''[Optional] If selected, the number of sequences per alignment
                    is not calculated after first trimming, and, thus, new deletions from the whole
                    dataset are prevented. However, your whole dataset might end up having
                    less alignments per sequence than proposed at `-w` option.''')
missingdata.add_argument('-l', '--stop_lib',
                    metavar = "",
                    type    = int,
                    choices = list(range(1, len(list( STOP_CODON_TABLE )) + 1 )),
                    default = 1,
                    help    = '[Optional] Stop library [Default = 1]')
missingdata.add_argument('-p','--print_stop_lib',
                    action="store_true",
                    help='[Optional] If selected, print stop libraries and exit')
missingdata.add_argument('-a','--add_deletion',
                    metavar = "",
                    default = None,
                    type    = str,
                    help    = '''[Optional] File in plain text with a list of 
                                species to delete on each alignment [default = None]''')
missingdata.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default = 1]')
missingdata.add_argument('-s','--suffix', 
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

pairedblast_required = pairedblast.add_argument_group('required arguments')
pairedblast_required.add_argument('-t','--taxonomy',
                    metavar="",
                    default = None,
                    required= True,
                    help='Taxonomy file. Format in csv: "[sequence name],[group]"')

pairedblast_positional = pairedblast.add_argument_group('positional arguments')
pairedblast_positional.add_argument('sequences', 
                    metavar='exons',
                    nargs="+",
                    type=str,
                    help='File names with sequences. If these are aligned, an unalignment process is performed')
pairedblast._action_groups.reverse()

# pairedblast ------------------------------------------------------

# bold search ------------------------------------------------------
boldsearch = subparsers.add_parser('bold',
                                    help = "Match sequences against the BOLD database",
                                    formatter_class = argparse.RawDescriptionHelpFormatter, 
                                    description="""


                Wrapper of both BOLD for species identification
                                    from COI sequences
- Host:
    BOLD: http://www.boldsystems.org/index.php/Ids_xml

Example:

    $ fishlifeqc bold [sequence] -t [taxonomy file]

        The taxnomy file is CSV-formated and must contain the following:

            names               spps
            [sequence header 1],[species names of header 1]
            [sequence header 2],[species names of header 2]
            [sequence header 3],[species names of header 3]
             ...                 ...                    

""")

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
# boldsearch.add_argument('-n','--ncbi',
#                         action='store_true',
#                         help=' If selected, BLASTn is used to identify species')
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

bold_required = boldsearch.add_argument_group('required arguments')
bold_required.add_argument('-t','--taxonomy',
                        metavar="file",
                        default = None,
                        required= True,
                        help='Taxonomy file. Format in csv: [sequence name],[species name]')

bold_positional = boldsearch.add_argument_group('positional arguments')
bold_positional.add_argument('sequence', 
                        metavar='COI file',
                        type=str,
                        help='''File name with the COI sequences. 
                                If these are aligned, an 
                                unalignment process is performed''')
boldsearch._action_groups.reverse()
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
tlike.add_argument('-f','--format',
                    metavar="",
                    type= str,
                    default= "newick",
                    help='[Optional] Tree format [Default: newick]') 
tlike.add_argument('-o','--outfile',
                    metavar="",
                    type= str,
                    default= "t_like_table.csv",
                    help='[Optional] Out filename [Default: t_like_table.csv]')
tlike.add_argument('-n', '--threads',
                    metavar = "",
                    type    = int,
                    default = 1,
                    help    = '[Optional] number of cpus [Default: 1]')   


tlike_tree_mod = tlike.add_argument_group('modification of input trees')
tlike_tree_mod.add_argument('-g','--outgroup',
                    metavar="",
                    type= str,
                    help='[Optional] Outgroup taxa [Default: None]')
tlike_tree_mod.add_argument('-l','--coll_bylen',
                    action="store_true",
                    help='''[Optional] If selected, collapse internal branches by length''')
tlike_tree_mod.add_argument('-u','--ucoll_bysupp',
                    action="store_false",
                    help='''[Optional] If selected, internal branches are not collapsed by support value''')
tlike_tree_mod.add_argument('-L', '--min_len',
                    metavar = "",
                    type    = float,
                    default = 0.000001,
                    help    = '''[Optional] minimun branch length to collapse 
                    internal branch. This value is also used to 
                    define T-like clades [Default: 0.000001]''')
tlike_tree_mod.add_argument('-S', '--min_supp',
                    metavar = "",
                    type    = float,
                    default = 0,
                    help    = '[Optional] minimun support value to collapse internal branch [Default: 0]')

tlike_required = tlike.add_argument_group('required arguments')
tlike_required.add_argument('-t','--taxonomyfile',
                    metavar="",
                    type= str,
                    help='Taxonomy file') 

tlike_positional = tlike.add_argument_group('positional arguments')
tlike_positional.add_argument('treefiles',
                    nargs="+",
                    help='Filenames')

tlike._action_groups.reverse()

# bl ------------------------------------------------------
bl = subparsers.add_parser('bl',
                            help = "Branch length ratios and correlations",
                            formatter_class = argparse.RawDescriptionHelpFormatter, 
                            description="""

                Branch length rations and pearson correlations

Both metrics are obtained by comparing a constrained
trees with a pruned reference tree.

Branch length ratios compares terminal branch lenghts. 
Pearson correlation uses branch lengths from root to tips.

Example:

    * Standard usage:

        $ fishlifeqc bl [sequences] -t [reference tree]

    * Codon aware running:

        $ fishlifeqc bl [sequences] -t [reference tree] -c 


""")

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

bl_raxml = bl.add_argument_group('RAxML constrained tree parameters')
bl_raxml.add_argument('-e','--evomol',
                 metavar="",
                 type= str,
                 default = 'GTRGAMMA',
                 help='[Optional] RAxML evol. model for constrained tree inference [Default: GTRGAMMA]')
bl_raxml.add_argument('-c','--codon_aware',
                 action="store_true",
                 help='[Optional] If selected, codon partition file is added')
bl_raxml.add_argument('-i', '--iterations',
                metavar = "",
                type    = int,
                default = 1,
                help    = '[Optional] Number of iterations for MLEs [Default: 1]')

bl_corr = bl.add_argument_group('correlation thresholds')
bl_corr.add_argument('-r', '--max_ratio',
                 metavar = "",
                 type    = float,
                 default = 5.0,
                 help    = '[Optional] Constrained/pruned tree maximum ratio allowed [Default: 5]')
bl_corr.add_argument('-p', '--min_pearson',
                 metavar = "",
                 type    = float,
                 default = 0.5,
                 help    = '[Optional] Minimum pearson correlation allowed [Default: 0.5]')

bl_required = bl.add_argument_group('required arguments')
bl_required.add_argument('-t','--reference',
                metavar="tree",
                type= str,
                default=None,
                required= True,
                help='Reference tree file in newick format [Default: None]')

bl_positional = bl.add_argument_group('positional arguments')
bl_positional.add_argument('sequences',
                nargs="+",
                type=str,
                help='Filenames')

bl._action_groups.reverse()
# bl ------------------------------------------------------


# para ------------------------------------------------------
para = subparsers.add_parser('para',
                            help = "Test paraphyly with AU tests",
                            formatter_class = argparse.RawDescriptionHelpFormatter, 
                            description="""

                AU test on paraphyletic groups

This command is designed to statiscially test (i.e., AU-test) the presence of a paraphyletic
group in a gene tree. If AU-test accept the alternative hypothesis (forced monophyly, see below),
it means the paraphyly statistically affect the topology of the whole gene tree.

Two hypothesis are tested: i) orignal tree (null) and ii) a constrained tree (alternative). 
The constrained tree is constructed by forcing the monophyly of a given paraphyletic
group on each gene tree. This given paraphyletic group is the most frequent group from 
the wholedataset. If a gene tree does not have this paraphyletic group, the second most 
frequent is used, and so on. 

Other options for selecting a paraphyletic group inside a gene tree are also available 
(see below options)

Example:

    * Standard usage:

        $ fishlifeqc para -A .fasta -T .tree  -t [taxonomy file]
        
            Where '-A' and '-T' indicate file extensions for alignments and 
            trees, correspondingly.

            Taxonomy file is CSV-formated and must contain the following
            struncture:

                # names       groups               
                [tip name 1],[group of tip name 1]
                [tip name 2],[group of tip name 2]
                [tip name 3],[group of tip name 3]
                ...          ...

    * Codon aware running:

        $ fishlifeqc para -t [taxonomy file] -c
""",
epilog= """
Note: Since constrained trees are reconstructed with RAxML internally, it is 
recommended to use RAxML trees as input trees to make better comparisons on 
site likelihood estimations (needed for AU tests with CONSEL)
""")
para.add_argument('-f','--forceall',
                   action="store_true",
                   help='''[Optional] If selected, force monophyly 
                            of all paraphyletic groups in a gene tree''')
para.add_argument('-g', '--test_group',
                metavar = "",
                type    = str,
                default = None,
                help    = '''[Optional] Group to force monophyly. If this group is not 
                            present in a given gene tree, this gene tree is skipped from 
                            analyses [Default: None]''')
para.add_argument('-o','--outfile',
                   metavar="",
                   type= str,
                   default= "au_tests.tsv",
                   help='[Optional] Out filename [Default: au_tests.tsv]')
para.add_argument('-n', '--threads',
                metavar = "",
                type    = int,
                default = 1,
                help    = '[Optional] number of cpus [Default: 1]')


para_raxml = para.add_argument_group('RAxML constrained tree parameters')
para_raxml.add_argument('-e','--evomol',
                 metavar="",
                 type= str,
                 default = 'GTRGAMMA',
                 help='[Optional] RAxML evol. model for constrained tree inference [Default: GTRGAMMA]')
para_raxml.add_argument('-c','--codon_aware',
                 action="store_true",
                 help='[Optional] If selected, codon partition file is added')
para_raxml.add_argument('-i', '--iterations',
                metavar = "",
                type    = int,
                default = 1,
                help    = '[Optional] Number of iterations for MLEs [Default: 1]')


para_tree_mod = para.add_argument_group('modification of input trees')
para_tree_mod.add_argument('-l','--coll_bylen',
                    action="store_true",
                    help='''[Optional] If selected, collapse internal branches by length''')
para_tree_mod.add_argument('-L', '--min_len',
                    metavar = "",
                    type    = float,
                    default = 0.000001,
                    help    = '''[Optional] minimun branch length to collapse 
                    internal branch [Default: 0.000001]''')
para_tree_mod.add_argument('-s','--coll_bysupp',
                    action="store_true",
                    help='''[Optional] If selected, internal branches are collapsed by support value''')
para_tree_mod.add_argument('-S', '--min_supp',
                    metavar = "",
                    type    = float,
                    default = 0,
                    help    = '[Optional] minimun support value to collapse internal branch [Default: 0]')

para_required = para.add_argument_group('input files')
para_required.add_argument('-p','--path',
                    metavar="",
                    type= str,
                    default=".",
                    help="Directory with trees and alignments [Default = '.']") 
para_required.add_argument('-A','--aln_ext',
                    metavar="",
                    type= str,
                    default=".fasta",
                    # required = True,
                    help="Alignment file extension [Default = '.fasta']") 
para_required.add_argument('-T','--tree_ext',
                    metavar="",
                    type= str,
                    default=".tree",
                    # required = True,
                    help="Tree file extension [Default = '.tree']") 
para_required.add_argument('-t','--taxonomyfile',
                    metavar="",
                    type= str,
                    required = True, 
                    help='Taxonomy file')

para._action_groups.reverse()
# para ------------------------------------------------------

# srh ------------------------------------------------------
srh = subparsers.add_parser('srh',
                            help = "Test stationarity, reversibility and homogeneity",
                            formatter_class = argparse.RawDescriptionHelpFormatter, 
                            description="""
        
    Assessing Stationarity, Reversibility and Homogeneity (SRH) assumptions

        Symmetry tests based on:
            * Bowker (1948) DOI: 10.1080/01621459.1948.10483284
            * Stuart (1955)  DOI: 10.2307/2333387
            * Ababneh et al. (2006) DOI: 10.1093/bioinformatics/btl064
        
        Code for maximum Match Pair Test of Symmetry based on:
            * Naser-Khdour et al. (2019) DOI: 10.1093/gbe/evz193

Examples:

    * Standard usage:

        $ %(prog)s [sequences]

    * Specify type of symmetry test:

        $ %(prog)s [sequences] -t sym

        note: there are three types: 'sym', 'mar', 'int'.
              'sym' (symmetry) is usually enough to test SRH assumptions.
              'mar' (marginal) approaches stationarity.
              'int' (internal) approaches homogeneity.

    * Codon aware: 

        $ %(prog)s [sequences] -c

        note: A partition file with position
              that passed test will be created
""")

srh.add_argument('sequences', 
                    metavar='exons',
                    nargs="*",
                    type=str,
                    help='''File names with sequences''')
srh.add_argument('-t', '--symtype',
                 metavar = "",
                 type    = str,
                 choices = ['sym', 'mar', 'int'],
                 default = 'sym',
                 help    = '''[Optional] symmetry test type. 
                             Option are 'sym', 'mar' and 'int' [Default = sym]''')
srh.add_argument('-p','--pval', 
                  metavar="",
                  type = float,
                  default = 0.05,
                  help='''[Optional] Threshold for P-values [Default: %s]''' % 0.05)
srh.add_argument('-r','--raxml',
                 action="store_false",
                 help='''[Optional] If selected, partitions are raxml-formated''')
srh.add_argument('-c','--codon_aware',
                  action="store_true",
                  help='[Optional] If selected, tests are done by codons')                 
srh.add_argument('-w','--write_bad',
                  action="store_true",
                  help='[Optional] If selected, not passing sequences are written')
srh.add_argument('-d','--trim_seqs',
                  action="store_true",
                  help='[Optional] If selected, trim not passing codon positions')
srh.add_argument('-a','--is_aa',
                  action="store_true",
                  help='[Optional] If selected, sequences are aminoacids')
srh.add_argument('-s','--suffix', 
                  metavar="",
                  default = "SymTest",
                  type= str,
                  help='[Optional] Suffix for outfile name [Default: %s]' % "SymTest")
srh.add_argument('-n', '--threads',
                  metavar = "",
                  type    = int,
                  default = 1,
                  help    = '[Optional] number of cpus [Default = 1]')
# srh ------------------------------------------------------


def main():

    wholeargs = parser.parse_args()

    if wholeargs.subcommand == "mdata":
        # print(wholeargs)
        init_class = Missingdata(
                fastas = wholeargs.sequences,
                htrim  = wholeargs.coverage, 
                vtrim  = wholeargs.edges, 
                itrim  = wholeargs.internal,
                min_sequences_per_aln = wholeargs.min_seqs_per_aln,
                min_alns_per_sequence = wholeargs.min_alns_per_seq,
                outputsuffix = wholeargs.suffix, 
                codon_aware  = wholeargs.codon_aware, # default false
                unadjusted= wholeargs.unadjusted,
                stop_opt = wholeargs.stop_lib,
                custom_deletion_list=wholeargs.add_deletion,
                threads  = wholeargs.threads,
            )

        if wholeargs.print_stop_lib:
            init_class.print_stop_lib()
            exit()

        if not wholeargs.sequences:
            print(missingdata.format_help())
            exit()

        init_class.run()

    elif wholeargs.subcommand == "rblast":
        from fishlifeqc.pairedblast import Pairedblast

        Pairedblast(
            sequences = wholeargs.sequences,
            taxonomy  = wholeargs.taxonomy,
            threads   = wholeargs.threads,
            outname   = wholeargs.out,
            threshold = wholeargs.identity
        ).run()

    elif wholeargs.subcommand == "bold":
        from fishlifeqc.boldsearch import Boldesearch
        
        Boldesearch(
            sequence           = wholeargs.sequence,
            bolddatabase       = wholeargs.bold_db,
            make_blast         = False, # wholeargs.ncbi, (API DEPRECATED)
            quiet              = wholeargs.quiet,
            taxonomyfile       = wholeargs.taxonomy,
            removeintermediate = wholeargs.keep,
            threshold          = wholeargs.threshold,
            outfile            = wholeargs.out
        ).id_engine()

    elif wholeargs.subcommand == "tlike":
        # print(wholeargs)
        from fishlifeqc.t_like import TreeExplore

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
        # print(wholeargs)
        from fishlifeqc.bl import BLCorrelations

        BLCorrelations(
            species_tree_file = wholeargs.reference,
            sequences         = wholeargs.sequences,
            ratio_threshold   = wholeargs.max_ratio,
            pearson_threshold = wholeargs.min_pearson,
            iterations        = wholeargs.iterations,
            evomodel          = wholeargs.evomol,
            prefix            = wholeargs.prefix,
            codon_partition   = wholeargs.codon_aware,
            threads           = wholeargs.threads
        ).BrLengths()

    elif wholeargs.subcommand == 'para':
        # print(wholeargs)
        from fishlifeqc.monophyly  import Monophyly

        Monophyly(
            path            = wholeargs.path,
            fasta_extension = wholeargs.aln_ext,
            tree_extension  = wholeargs.tree_ext,
            taxonomyfile    = wholeargs.taxonomyfile,
            
            collapsebylen   = wholeargs.coll_bylen,   # for collapse
            minlen          = wholeargs.min_len,      # for collapse
            collpasebysupp  = wholeargs.coll_bysupp, # for collapse
            minsupp         = wholeargs.min_supp,     # for collapse
            force_all       = wholeargs.forceall,   # if True, it will force mono to ALL para
            tgroup          = wholeargs.test_group, # target group
            evomodel        = wholeargs.evomol,
            codon_partition = wholeargs.codon_aware,
            iterations      = wholeargs.iterations,
            out_report      = wholeargs.outfile,
            threads         = wholeargs.threads
        ).run()

    elif wholeargs.subcommand == 'srh':
        from fishlifeqc.symtests import SymTests # expensive import

        SymTests(
            sequences   = wholeargs.sequences,
            codon_aware = wholeargs.codon_aware,
            suffix      = wholeargs.suffix,
            symtype     = wholeargs.symtype,
            nexusformat = wholeargs.raxml,
            pval        = wholeargs.pval,
            threads     = wholeargs.threads,
            write_bad   = wholeargs.write_bad,
            isaminoacid = wholeargs.is_aa,
            trim_seqs   = wholeargs.trim_seqs,
        ).main()

if __name__ == "__main__":
    main()
