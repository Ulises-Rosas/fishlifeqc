#!/usr/bin/env python3

import os
import argparse
from multiprocessing import Pool

import fishlifeseq
from fishlifeqc.utils import fas_to_dic

class Concatenate:
    def __init__(self,
                 alignments = None,
                 supermatrixname = None,
                 partitionsname = None,
                 iscodon = False,
                 nexusformat = True,
                 threads = 3):

        self.alignments = alignments
        self.threads    = threads
        self.supermatrixname = supermatrixname
        self.partitionsname = partitionsname
        self.iscodon = iscodon
        self.nexusformat = nexusformat

        # placeholder
        self.uniqueheaders = []

    def completeseq(self, myfile):

        fasta = fas_to_dic(myfile)

        seqheaders = list(fasta.keys())
        alnlen     = len(fasta[seqheaders[0]])
        leftseqs   = self.uniqueheaders - set(seqheaders) 

        if leftseqs:
            for ls in leftseqs:
                fasta[ls] = '-'*alnlen

        return (myfile, fasta, alnlen)

    def formating(self, exon, pos1, pos2):

        exonb = os.path.basename(exon)
        out = []
        if self.iscodon:
            if self.nexusformat:
                # "  charset %s_pos1 = %s: %s-%s\\3;\n" % (exonb, exon, pos1    , pos2),
                out.extend([
                    "  charset %s_pos1 = %s-%s\\3;\n" % (exonb, pos1    , pos2),
                    "  charset %s_pos2 = %s-%s\\3;\n" % (exonb, pos1 + 1, pos2),
                    "  charset %s_pos3 = %s-%s\\3;\n" % (exonb, pos1 + 2, pos2)
                    ])
            else:
                out.extend([
                    "DNA, %s_pos1 = %s-%s\\3\n" %(exonb, pos1    , pos2),
                    "DNA, %s_pos2 = %s-%s\\3\n" %(exonb, pos1 + 1, pos2),
                    "DNA, %s_pos3 = %s-%s\\3\n" %(exonb, pos1 + 2, pos2)
                    ])

        else:
            if self.nexusformat:
                out.append(
                    "\tcharset %s = %s-%s;\n" % (exonb, pos1, pos2)
                )
            else:
                out.append(
                   "DNA, %s = %s-%s\n" %(exonb, pos1, pos2)
                )

        return out

    def reducedict(self, fastainfo):

        # fastainfo = extendedfastas
        merged  = {}
        pos  = 1
        partitions = []

        for exon,fasta_dict,itslen in fastainfo:

            partitions.extend( 
                self.formating(exon, pos, pos + itslen - 1)
            )
            pos  += itslen

            for k,v in fasta_dict.items():
                if not merged.__contains__(k):

                    merged[k] = v
                else:
                    merged[k] += v


        with open(self.supermatrixname, 'w') as f:
            for k,v in merged.items():
                f.write( "%s\n%s\n" % (k,v) )

        with open(self.partitionsname, 'w') as f:

            if self.nexusformat:
                f.write('#nexus\n')
                f.write('begin sets;\n')

            f.writelines(partitions)

            if self.nexusformat:
                f.write('end;\n')

    def run(self):

        allheaders  = []
        with Pool(processes = self.threads) as p:

            preallheaders = [*p.map(fishlifeseq.headers, self.alignments)]

            for i in preallheaders:
                allheaders.extend(i)

            self.uniqueheaders = set(allheaders)
            extendedfastas     = [*p.map(self.completeseq, self.alignments)]

            self.reducedict(extendedfastas)


# self = Concatenate(alignments=['../simple1.txt', '../simple2.txt'], iscodon=True, nexusformat=False)


def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

                        Concatenate sequences

    * Concatenate without creating a codon partition file (default):

        $ concatenate.py [exon files]

    * Concatenate and create a codon partition file:

        $ concatenate.py [exon files] --codon_aware

    * Let the partition file take a raxml-based format:

        $ concatenate.py [exon files] --raxml
""")
    parser.add_argument('filenames',
                        metavar = 'exons',
                        nargs="+",
                        help='Filenames with sequences')
    parser.add_argument('-c','--codon_aware',
                        action="store_true",
                        help='''[Optional] If selected, codon partitions file is created''')
    parser.add_argument('-r','--raxml',
                        action="store_false",
                        help='''[Optional] If selected, partition file is raxml-formated''')
    parser.add_argument('-o','--out_name',
                        metavar="",
                        type= str,
                        default= "mysupermatrix.txt",
                        help='[Optional] Concatenate output name [Default = "mysupermatrix.txt"]') 
    parser.add_argument('-p','--partition_name',
                        metavar="",
                        type= str,
                        default= "mypartitions.txt",
                        help='[Optional] Name of partition file [Default = "mypartitions.txt"]') 
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')                        
    args = parser.parse_args()
    return args

def main():
    args = getOpts()
    # print(args)

    Concatenate(
        alignments      = args.filenames ,
        supermatrixname = args.out_name,
        partitionsname  = args.partition_name,
        nexusformat     = args.raxml, # default: True
        iscodon         = args.codon_aware, # default: False
        threads         = args.threads
        ).run()

if __name__ == "__main__":
    main()
