import sys
import argparse
from multiprocessing import Pool
from fishlifeqc.utils import isfasta
from boldminer.utils import fas_to_dic

FAILEDTOTRIM = "failed_to_trim.txt"

class Missingdata:

    def __init__(self,
                fastas, 
                htrim = 0.6, 
                vtrim = 0.5, 
                outputsuffix = "_trimmed", 
                trimedges = False,
                threads = 1):

        self.fasta = fastas
        self.htrim = htrim
        self.vtrim = vtrim
        self.outputsuffix = outputsuffix
        self.removeedges  = trimedges
        self.threads      = threads

    def check_aln(self, aln, filename):
        """
        it gets dictionary 
        data structure
        """
        # aln = alignment
        lengths = set([len(v) for _,v in aln.items()])

        if lengths.__len__() > 1:
            sys.stderr.write("'%s' file has sequences with different lengths\n" % filename)
            sys.stderr.flush()
            return None

        else:
            return lengths.pop()

    def makebarplots(self, simplelist):
        # vals = lftiout
        import matplotlib.pyplot as plt

        vals = simplelist
        breaks = [n for n,_  in enumerate(vals) ]
        
        plt.figure(figsize=(10,5))
        plt.bar(breaks, vals,  align='center')
        plt.xlabel('Positions')
        plt.title('Gap percentage')
        plt.show()
        plt.close()

    def sequencecompleteness(self, aln, seqlength):
        # remove those
        # sequences when
        # the gap perce is more
        # than the threshold value
        out = {}
        for k,v in aln.items():
            if v.count('-')/seqlength <= self.htrim:
                out.update( {k:v}  )
        return out

    def trimedges(self, aln, seqlength):
        # remove those
        # columns when
        # the gap perce is more
        # than the threshold value

        # threshold = 0.5
        headerlength = aln.__len__()
        
        lfti = 0

        while lfti < seqlength:
            gap_perc = [v[lfti] for _,v in aln.items()].count('-')/headerlength

            if gap_perc <= self.vtrim:
                break
            lfti += 1

        rgti = seqlength - 1

        while rgti >= 0:
            gap_perc = [ v[rgti] for _,v in aln.items() ].count('-')/headerlength

            if gap_perc <= self.vtrim:
                break
            rgti -= 1

        if rgti > lfti:
            return { k:v[ lfti:rgti + 1 ] for k,v in aln.items() }
        else:
            return None

    def writeresults(self, obj, name):

        outname = name + self.outputsuffix

        with open(outname, 'w') as f:

            for k,v in obj.items():
                f.write( "%s\n%s\n" % (k,v))

    def trimiterator(self, fasta):

        if not isfasta(fasta):
            return None

        alignment = fas_to_dic(fasta)
        seqlength = self.check_aln(alignment, fasta)

        if not seqlength:
            return None

        # horizontal trim
        trimmed = self.sequencecompleteness(
                    aln       = alignment,
                    seqlength = seqlength
                )

        if self.removeedges:
            # vertical trim
            trimmed = self.trimedges(
                        aln       = trimmed, 
                        seqlength = seqlength
                    )

        if trimmed:
            self.writeresults(trimmed, fasta)
            return fasta
        else:
            return None

    def run(self):
        with Pool(processes = self.threads) as p:
            passed = [*p.map(self.trimiterator, self.fasta)]

        failed = set(self.fasta) - set(list(filter(None, passed)))

        if failed:
            with open(FAILEDTOTRIM, 'w') as f:
                for i in failed:
                    f.write( "%s\n" % i)
