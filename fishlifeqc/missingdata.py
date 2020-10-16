import sys
import argparse
from multiprocessing import Pool
from fishlifeqc.utils import isfasta
from fishlifeqc.utils import fas_to_dic

FAILEDTOTRIM = "failed_to_trim.txt"

class Missingdata:

    def __init__(self,
                fastas = None, 
                htrim = 0.6, 
                vtrim = 0.5, 
                outputsuffix = "_trimmed", 
                # trimedges   = False,
                horizontal  = True,
                codon_aware = False,
                mtlib       = False,
                quiet       = False,
                threads     = 1):

        self.fasta = fastas
        self.htrim = htrim
        self.vtrim = vtrim
        self.outputsuffix = outputsuffix
        self.horizontal   = horizontal
        # self.removeedges  = trimedges
        #   |
        #    -> if removeedges:
        self.codon_aware  = codon_aware
        #   |
        #    -> if codon_aware:
        self.mtlib        = mtlib
        self.threads      = threads

        # stops for mitchondrial DNA
        self.mt_stop_codons=["TAA", "TAG", "AGA", "AGG"]
        # stops for standard genetic code
        self.std_stop_codons=["TAA", "TAG", "TGA"]

        self.quiet = quiet
        # placeholder
        # self.count = 0

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

    def sequencecompleteness(self, aln):
        # remove those
        # sequences when
        # the gap perce is more
        # than the threshold value
        seqlength = len(aln[ list(aln)[0] ])
        out = {}
        for k,v in aln.items():
            if v.count('-')/seqlength <= self.htrim:
                out.update( {k:v}  )
        return out

    def trimedges(self, aln):
        # remove those
        # columns when
        # the gap perce is more
        # than the threshold value

        # threshold = 0.5
        seqlength = len(aln[ list(aln)[0] ])
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

    def is_stop(self, codon):

        if self.mtlib:
            return codon in self.mt_stop_codons
        
        else:
            return codon in self.std_stop_codons

    def filter_stopcodons(self, aln):
        """
        stop codon checker
        """
        
        seqlength = len(aln[ list(aln)[0] ])

        out = {}
        for k,v in aln.items():
            stopcount = 0
            mystr     = ""
            
            for c in range(0, seqlength, 3):                
                codon = v[c:c+3]

                if self.is_stop(codon):
                    stopcount += 1
                    mystr     += "NNN"

                else:
                    mystr += codon

                if stopcount > 1:
                    break
                    
            if stopcount <= 1:
                out[k] = mystr

        return out
        
    def trimedges_codonaware(self, aln):
        
        headerlength = aln.__len__()
        seqlength = len(aln[ list(aln)[0] ])
        
        lfti = 0
        for i in range(0, seqlength, 3):

            F = [v[i]     for _,v in aln.items()].count('-')
            S = [v[i + 1] for _,v in aln.items()].count('-')
            T = [v[i + 2] for _,v in aln.items()].count('-')
            perc = (F + S + T)/(headerlength*3)

            if perc <= self.vtrim:
                lfti = i
                break
                
        rgti = 0
        for i in reversed(range(0, seqlength, 3)):

            F = [v[i]     for _,v in aln.items()].count('-')
            S = [v[i + 1] for _,v in aln.items()].count('-')
            T = [v[i + 2] for _,v in aln.items()].count('-')

            perc = (F + S + T)/(headerlength*3)

            if perc <= self.vtrim:
                rgti = i
                break
                
        if rgti > lfti:
            trimmed = { k:v[ lfti:rgti ] for k,v in aln.items() }
             # horizontal trimming (check stop codons)
            return self.filter_stopcodons(trimmed)
        else:
            return None
            
    def close_gaps(self, aln: dict, is_codon_aware: bool = True) -> dict:
        """
        It closes the empty columns formed by 
        `self.sequencecompleteness` function
        """
        # aln = {'a':'---cat--a', 'b':'---cat--a'}
        # is_codon_aware = not True

        offset = 3 if is_codon_aware else 1
        seqlength = len(aln[ list(aln)[0] ])

        idxs = []
        for c in range(0, seqlength, offset):
            mystr = ""

            for k,v in aln.items():
                mystr += v[c:c+offset]

            is_allgap = 1 == mystr.count('-')/len(mystr)

            if not is_allgap:
                idxs.append(c)

        out = {}
        for k,v in aln.items():
            mystr  = ""

            for c in range(0, seqlength, offset):
                if c in idxs:
                    mystr += v[c:c+offset]

            out[k] = mystr
            
        return out   
    
    def writeresults(self, obj, name):

        outname = name + self.outputsuffix

        with open(outname, 'w') as f:

            for k,v in obj.items():
                f.write( "%s\n%s\n" % (k,v))

    def trimiterator(self, fasta):

        # fasta = self.fasta[0]
        if not isfasta(fasta):
            return None

        alignment = fas_to_dic(fasta)
        seqlength = self.check_aln(alignment, fasta)

        if not seqlength:
            return None

        # TODO: delete extreme cases in sequence length 

        if not self.codon_aware:
            # vertical trim, length change
            trimmed = self.trimedges(aln = alignment)
        else:
            # vertical trim, length change &
            # horizontal trimming (check stop codons)
            trimmed = self.trimedges_codonaware(aln = alignment)

        if self.horizontal:
            # horizontal trim, not length change
            trimmed = self.sequencecompleteness(aln = trimmed)
        
        if trimmed:
            # close the gaps
            trimmed = self.close_gaps(trimmed, self.codon_aware)

        if not self.quiet:
            # self.count += 1
            # totallength = len(self.fasta)
            # # glob_count += 1
            # print(self.count)
            sys.stderr.write("Processed: %s\n" % fasta )
            sys.stderr.flush()

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
