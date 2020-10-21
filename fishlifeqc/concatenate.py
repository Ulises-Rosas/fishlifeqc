
import os
import argparse
from multiprocessing import Pool

import fishlifeseq
from fishlifeqc.utils import fas_to_dic

class Concatenate:
    def __init__(self,
                 alignments,
                 supermatrixname,
                 partitionsname,
                 threads = 3):

        self.alignments = alignments
        self.threads    = threads
        self.supermatrixname = supermatrixname
        self.partitionsname = partitionsname

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

    def reducedict(self, fastainfo):

        merged  = {}
        pos  = 1
        ROWPART = "DNA, %s=%s-%s\n"
        partitions = []

        for exon,fasta_dict,itslen in fastainfo:

            exonb = os.path.basename(exon)
            partitions.append( 
                ROWPART % (exonb, pos, pos + itslen - 1) 
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
            f.writelines(partitions)

    def run(self):

        allheaders  = []
        with Pool(processes = self.threads) as p:

            preallheaders = [*p.map(fishlifeseq.headers, self.alignments)]

            for i in preallheaders:
                allheaders.extend(i)

            self.uniqueheaders = set(allheaders)
            extendedfastas     = [*p.map(self.completeseq, self.alignments)]

            self.reducedict(extendedfastas)

def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

                    Concatenate sequences

    * Example:

        $ concatenate.py [exon files]
""")
    parser.add_argument('filenames',
                        metavar = 'exons',
                        nargs="+",
                        help='Filenames with sequences')
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

    Concatenate(
        alignments = args.filenames ,
        supermatrixname  = args.out_name,
        partitionsname = args.partition_name,
        threads = args.threads
        ).run()

if __name__ == "__main__":
    main()

# myfiles = ['data/COI.NT_aligned.fasta',
#            'data/E0537.NT_aligned.fasta',
#            'data/E1718.NT_aligned.fasta'
#            ]

# Concatenate( 
#     alignments = myfiles,
#     supermatrixname= "mysupermatrix.txt",
#     partitionsname = "mypartitions.txt",
#     threads    = 3
# ).run()

