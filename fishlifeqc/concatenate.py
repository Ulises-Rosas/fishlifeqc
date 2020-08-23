
import os
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

