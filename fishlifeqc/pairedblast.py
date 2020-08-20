
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
from fishlifeqc.utils import isfasta

myos = sys.platform

if myos == 'darwin':
    BLASTN      = 'blastn_Darwin_64bit'
    MAKEBLASTDB = 'makeblastdb_Darwin_64bit'

elif myos == 'linux' or  myos == 'linux2':
    BLASTN      = 'blastn_Linux_64bit'
    MAKEBLASTDB = 'makeblastdb_Linux_64bit'

elif myos == 'linux' or  myos == 'linux2':
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
                threads   = 1,
                threshold = 95):
        """
        makeblastdb –in mydb.fsa –dbtype nucl
        blastn -db db -query query -out  out -outfmt 5 -perc_identity 95
        """
        self.threads   = threads
        self.sequences = sequences
        self.taxonomy  = taxonomy
        self.threshold = threshold

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
                seq, group, _ = l.strip().split(",")
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

    def iteratedbproduction(self):

        with Pool(processes = self.threads) as p:
            mydbs = [*p.map(self.generatedb, self.sequences)]            
        # mydbs = [*map(self.generatedb, 
        # ["../data/test_hyphen.txt", "../data/exon1.txt"])]

        try:
            os.remove(self.makeblastlogfil)
            os.remove(self.makeblastlogfil.split(".")[0] + ".perf")
        except FileNotFoundError:
            pass

        passed = list(filter(None, mydbs))
        failed = list(set(self.sequences) - set(passed))

        return (passed, failed)

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

            othermatch = self.checkoutgroups(q, tmpout, filename)

            if othermatch:
                out.extend(othermatch)

            os.remove(q)
            os.remove(tmpout)

        if not failedtoblast:
            [*map(os.remove, glob.glob(nohyphen+".n*"))]
            os.remove(nohyphen)

        else:
            for line in failedtoblast:
                sys.stderr.write("Error runnig: %s\n" % line)
                sys.stderr.flush()

        if not out:
            return None

        return list(set(out))

    def iterateblastn(self, passed):

        with Pool(processes = self.threads) as p:
            preout = [*p.map(self.blastn, passed)]

        # preout = [*map(self.blastn, 
        # ["../data/test_hyphen.txt", "../data/exon1.txt"])]

        preout = filter(None, preout)

        if preout:
            out  = []
            for i in preout:
                out += i

            return out

        else:
            return None
