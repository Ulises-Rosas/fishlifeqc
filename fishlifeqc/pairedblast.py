
import re
import os
import sys
import glob
import platform
import argparse
from multiprocessing import Pool
import xml.etree.ElementTree as ET


import runshell
import fishlifeseq
from fishlifeqc.utils import isfasta,fas_to_dic
from fishlifeqc.missingdata import Missingdata

myos = sys.platform

if myos == 'darwin':
    BLASTN      = 'blastn_Darwin_64bit'
    MAKEBLASTDB = 'makeblastdb_Darwin_64bit'

elif myos == 'linux' or  myos == 'linux2':
    BLASTN      = 'blastn_Linux_64bit'
    MAKEBLASTDB = 'makeblastdb_Linux_64bit'

elif myos == 'win32':
    BLASTN      = 'blastn.exe'
    MAKEBLASTDB = 'makeblastdb.exe'


def getqueries(prefix, filename):

    ourpattern = filename + "_query[0-9]+$"
    out = []

    if not prefix:
        prefix = "."

    for i in os.listdir(prefix):
        if re.findall(ourpattern, i):
            out.append(
                os.path.join(prefix, i)
                )

    return out


class Pairedblast:

    def __init__(self,
                sequences = None,
                taxonomy  = None,
                codon_aware = False,
                threads   = 1,
                outname = "mismatch_pairedblastn.txt",
                threshold = 95):
        """
        makeblastdb –in mydb.fsa –dbtype nucl
        blastn -db db -query query -out  out -outfmt 5 -perc_identity 95
        """
        self.threads      = threads
        self.sequences    = sequences
        self.taxonomy     = taxonomy
        self.codon_aware  = codon_aware
        self.threshold    = threshold
        self.outname      = outname
        self.blastfailure = outname + "_failed_to_makeblastdb"

        self.suffix       = "_rblastd"
        self.failedelete  = outname + "_failed_to_delete"

        self.blastnexec = BLASTN
        self.PREFIX     = "noaln_"
        self.blastline  = "{blastn} -db {db}  -query {query} -out {outp} -outfmt 5 -perc_identity {threshold}"
        
        self.makeblastdbexec = MAKEBLASTDB
        self.makeblastlogfil = "fishlife_makeblast.log"
        self.makeblastdbline = "{makeblastdb} -in {inp}  -dbtype nucl -logfile {log}" 

        self.ADDRESS = "BlastOutput_iterations/Iteration/Iteration_hits"

    @property
    def readtax(self):
        # self.taxonomy = "../data/taxonomy.txt"
        if self.taxonomy is None:
            return None

        out = {}
        with open(self.taxonomy, 'r') as f:
            for l in f.readlines():
                # print(l.strip().split(","))
                seq, group = l.strip().split(",")
                out[seq] = group
        
        return out

    def checktaxon(self, sequence):

        headers  = []

        with open(sequence, 'r') as f:
            filelines = f.readlines()
            for fl in filelines:
                if re.findall("^>", fl):
                    headers.append( fl.strip().replace(">", "") )

        headersuniq = list(set(headers))

        out = []
        for h in headersuniq:
            if not self.readtax.__contains__(h):
                out.append(h)
        # sequences not
        # represented at the
        # taxonomy file
        return out

    def checktaxonomyfile(self):
        # test = ["../data/test_hyphen.txt", "../data/exon1.txt"]
        # preout = [*map(self.checktaxon, test)]

        with Pool(processes = self.threads) as p:
            preout = [*p.map(self.checktaxon, self.sequences)]
        
        out = []
        for i in preout:
            out += i

        if out:
            sys.stderr.write("Check the following at the taxonomy file:\n\n")
            sys.stderr.flush()

            for i in set(out):
                sys.stderr.write("- %s\n" % i )
                sys.stderr.flush()
            exit()

    def parseblastnanes(self, blastnout):

        tree  = ET.parse(blastnout).getroot()

        outnames = []
        for child in tree.find(self.ADDRESS):
            outnames.append(child.find('Hit_def').text)

        return outnames

    def checkoutgroups(self, query, blastnout, filename):

        headquery = ""

        with open(query, 'r') as myfile:
            for i in myfile:
                if re.findall("^>", i):
                    headquery = i.strip().replace(">", "")
                    break
        
        targetgroup = self.readtax[headquery]
        outgroup    = []

        parsed = self.parseblastnanes(blastnout)

        for spps in parsed:
            tmpgroup = self.readtax[spps]

            if tmpgroup != targetgroup:
                outgroup += [(filename, spps, tmpgroup)]

        return outgroup

    def generatedb(self, sequence):

        if not isfasta(sequence):
            return None

        filename = os.path.basename(sequence)
        prefix   = os.path.dirname(sequence)
        nohyphen = os.path.join(prefix, self.PREFIX + filename)

        try:
            fishlifeseq.rm_hyphens(sequence, nohyphen)
        except ValueError:
            return None
        
        status = runshell.get(
                    self.makeblastdbline.format(
                        makeblastdb = self.makeblastdbexec,
                        inp         = nohyphen,
                        log         = self.makeblastlogfil
                        )
                    )

        if not status:
            return sequence

        else:
            return None

    def blastn(self, sequence):
        # sequence = "../data/test_hyphen.txt"

        filename = os.path.basename(sequence)
        prefix   = os.path.dirname(sequence)
        nohyphen = os.path.join(prefix, self.PREFIX + filename)
    
        status = fishlifeseq.splitfasta(nohyphen)

        if status:
            sys.stderr.write("Unable to split fastas at '%s' file\n" % filename)
            sys.stderr.flush()
            exit()

        queries = getqueries(prefix, filename)

        out = []
        failedtoblast = []

        for q in queries:
            tmpout = q + "_blastn"
            line   = self.blastline.format(
                            query     = q,
                            outp      = tmpout,
                            db        = nohyphen, 
                            threshold = self.threshold,
                            blastn    = self.blastnexec
                        )

            status = runshell.get(line)

            if status:
                failedtoblast.append(line)
                continue

            othermatch = self.checkoutgroups(q, tmpout, sequence)

            if othermatch:
                out.extend(othermatch)

            os.remove(q)
            os.remove(tmpout)

        if not failedtoblast:
            [*map(os.remove, glob.glob(nohyphen+".n*"))]
            os.remove(nohyphen)
            sys.stderr.write("Processed: %s\n" % sequence )
            sys.stderr.flush()        

        else:
            for line in failedtoblast:
                sys.stderr.write("Error runnig: %s\n" % line)
                sys.stderr.flush()

        if not out:
            return None

        return list(set(out))

    def constraint(self, myoutliers):

        myfile, badheaders = myoutliers

        myfasta  = fas_to_dic(myfile)

        newfasta = {}
        for h,s in myfasta.items():

            if not h.replace(">", "") in set(badheaders):
                newfasta[h] = s

        if not newfasta:
            sys.stderr.write("'%s' file got all mismatched\n" % myfile)
            sys.stderr.flush()
            return myfile

        with open( myfile + self.suffix, 'w') as f:
            # Close gaps :
            newfasta = Missingdata.close_gaps(newfasta, self.codon_aware)

            for h,s in newfasta.items():
                f.write("%s\n%s\n" % (h,s))
        return None

    def rebranding(self, myfile):
        myfasta  = fas_to_dic(myfile)

        with open( myfile + self.suffix, 'w') as f:
            for h,s in myfasta.items():
                f.write("%s\n%s\n" % (h,s))

    def run(self):

        self.checktaxonomyfile()

        with Pool(processes = self.threads) as p:

            ok_sequences = set(self.sequences)

            mydbs = [*p.map(self.generatedb, ok_sequences)]            

            try:
                os.remove(self.makeblastlogfil)
                os.remove(self.makeblastlogfil.split(".")[0] + ".perf")
            except FileNotFoundError:
                pass

            passed = list(filter(None, mydbs))
            failed = list(ok_sequences - set(passed))

            # return (passed, failed)
            # passed, failed = self.iteratedbproduction()

            if failed:
                ok_sequences -= set(failed)

                with open(self.blastfailure, "w") as f:
                    for i in failed:
                        f.write(i + "\n")

            if passed:

                preout = [*p.map(self.blastn, passed)]
                preout = list(filter(None, preout))

                # mismatches
                if preout:

                    outliers  = []
                    for po in preout:
                        outliers.extend(po)

                    badseqs = {}
                    with open(self.outname, 'w') as f:
                        f.write("exon,sample,group\n")

                        for exon,spps,group in outliers:

                            ok_sequences -= set([exon])
                            
                            if badseqs.__contains__(exon):
                               badseqs[exon] += [spps]
                            else:
                                badseqs[exon] = [spps]

                            f.write("%s,%s,%s\n" %  (exon, spps, group))

                    failed = [*p.map(self.constraint, tuple(badseqs.items()))]
                    failed = list(filter(None, failed))

                    if failed:
                        with open(self.failedelete, 'w') as f:
                            for i in failed:
                                f.write( i + "\n")

            [*p.map(self.rebranding, ok_sequences)]


# taxonomy  =  'data/taxonomy_genera.txt'
# sequences = ['data/mock1.txt',
#              'data/mock2.txt',
#              'data/mock3.txt']

# Pairedblast(sequences = sequences,
#             taxonomy  = taxonomy,
#             threads   = 4,
#             outname     = "mismatch_pairedblastn.txt",
#             threshold = 90).run()

