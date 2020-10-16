import os
import re
import argparse
from multiprocessing import Pool
from fishlifeqc.utils import fas_to_dic


def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

            Merge extracted information from one directory to another

    * Usage:
        $ merge.py [exon files] -m [map file] -n [ncpus]
        
        Note: Map file has the following structure:
                [locus1],[file name with locus1]
                [locus2],[file name with locus2]
                [locus3],[file name with locus3]

        Note 2: Exon files are have the this structure: "[locusX].*"
""")
    parser.add_argument('filenames',
                        nargs="+",
                        help='exons names')
    parser.add_argument('-m','--mapfile',
                        metavar="",
                        required=True,
                        type= str,
                        default= "mapfile.txt",
                        help='[Optional] Map file with locus names and file names')                    
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')

    args = parser.parse_args()
    return args

if __name__ == "__main__":

    args = getOpts()

    mapfile   = args.mapfile
    exonfiles = args.filenames

    mapfiles = []

    with open(mapfile, 'r') as f:
        for i in f.readlines():
            locus,tmpfile = i.strip().split(',')
            mapfiles.append((locus, tmpfile))

    def mergefasta(locus_tmpf):

        locus,tmpf = locus_tmpf
        locpatt    = re.compile("^%s\\." % locus)

        # matched_locus = []
        matched = ""

        for tolocus in exonfiles:
            
            exonfilebase   = os.path.basename(tolocus)
            if locpatt.match(exonfilebase):

                outfile        = os.path.join(tolocus + "_lily")                
                fasta_mainfile = fas_to_dic(file = tolocus)

                with open(outfile, 'w') as f:
                    for k1,v1 in fasta_mainfile.items():
                        f.write( "%s\n%s\n" % (k1,v1) )

                    for k2,v2 in fas_to_dic(file = tmpf).items():
                        f.write( "%s\n%s\n" % (k2,v2) )

                # print(tolocus)
                matched += tolocus
                break

        return matched

    def rebranding(restlocus):

        outfile = os.path.join(restlocus + "_lily")
        fastad  = fas_to_dic(file = restlocus)

        with open(outfile, 'w') as f:
            for k1,v1 in fastad.items():
                f.write( "%s\n%s\n" % (k1,v1) )

    with Pool(processes = args.threads) as p:
        # print(mapfiles)

        matched  = [*p.map(mergefasta, mapfiles)]
        restloci = set(exonfiles) - set(list(filter(None, matched)))
        # print(restloci)
        [*p.map(rebranding, restloci)]

