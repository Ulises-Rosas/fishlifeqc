import os
import re
import sys

import fishlifeseq
from fishlifeqc.utils import fas_to_dic
from boldminer.id_engine import id_engine 

class Taxonomychecks:
    
    def __init__(self,
                 taxonomyfile = None,
                 sequence = None,
                 specieslevel = True,
                 threads = 1):
        
        self.taxonomyfile = taxonomyfile
        self.sequence = sequence
        self.specieslevel = specieslevel
        self.threads = threads
        
    @property
    def sppstax(self):
        
        if self.taxonomyfile is None:
            sys.stderr.write('Please introduce a proper taxonomy file\n')
            sys.stderr.flush()
            exit()
            
        out = {}
        with open(self.taxonomyfile, 'r') as f:
            for l in f.readlines():
                # print(l.strip().split(","))
                seq, spps = l.strip().split(",")
                out[seq] = spps

        return out
    
    @property
    def getheaders(self):
        
        headers  = []
        with open(self.sequence, 'r') as f:
            filelines = f.readlines()
            for fl in filelines:
                if re.findall("^>", fl):
                    headers.append( fl.strip().replace(">", "") )

        return list(set(headers))
    
    def checktaxon(self):
        
        headersuniq = self.getheaders
                
        if self.specieslevel:
            readtaxonomy = self.sppstax
                    
        out = []
        for h in headersuniq:
            if not readtaxonomy.__contains__(h):
                out.append(h)
        # sequences not
        # represented at the
        # taxonomy file
        return out
    
        
    def check_spps_name(self):
        spps_pat  = "^[A-Z][a-z]+ [a-z]+$"
        anti_pat1 = "^[A-Z][a-z]+ sp[p]{0,1}$"
        
        headersuniq = self.getheaders
        
        myspps = {}
        notspps = []
        NAspps = []
        
        for i in headersuniq:
            tmp_spps = self.sppstax[i]
            
            if re.findall(spps_pat , tmp_spps) and\
               not re.findall(anti_pat1, tmp_spps):
                myspps[i] = tmp_spps
                
            else:
                if not tmp_spps or tmp_spps == "NA":
                    NAspps.append(i)
                    
                else:
                    notspps.append(tmp_spps)
                    
        return (myspps, notspps, NAspps)

class Boldesearch:

    def __init__(self,
                 sequence = None,
                 bolddatabase = 'COX1_SPECIES_PUBLIC',
                 make_blast = False,
                 quiet = False,
                 taxonomyfile = None,
                 removeintermediate = True,
                 threshold = 0.98,
                 outfile = None):

        self.sequence = sequence
        self.bolddatabase = bolddatabase
        self.make_blast = make_blast
        self.quiet    = quiet
        self.taxonomyfile = taxonomyfile
        self.threshold = threshold
        self.removeintermediate = removeintermediate
        self.outfile = outfile
        
    @property
    def dealignment(self):
        nohyphen = self.sequence + "_noaln"
        status = fishlifeseq.rm_hyphens(self.sequence, nohyphen)
        
        if status:
            sys.stderr.write("There were issues removing '-' chars at '%s' file\n" % self.sequence)
            sys.stderr.flush()
            exit()
        return nohyphen
    
    def checkspps(self):
        mydealn   = self.dealignment
        taxocheck = Taxonomychecks(taxonomyfile = self.taxonomyfile, sequence = mydealn)
    
        missingtaxons = taxocheck.checktaxon()
        
        if missingtaxons:
            sys.stderr.write("These sequences are not covered in '%s' file:\n\n" % self.taxonomyfile)            
            for i in missingtaxons:
                sys.stderr.write(" - %s\n" % i)
            sys.stderr.flush()            
            exit()
            
        myspps, notspps, NAspps = taxocheck.check_spps_name()
        
        if not self.quiet:
            
            if notspps:
                sys.stderr.write("These names are not binomial species names:\n\n")
                for i in notspps:
                    sys.stderr.write(" - %s\n" % i)

            if NAspps:
                sys.stderr.write("\nFurthermore, these sequences do not have species names:\n\n")
                for i in NAspps:
                    sys.stderr.write(" - %s\n" % i)

            if notspps or NAspps:
                sys.stderr.write("\nCheck your '%s' file\n" % self.taxonomyfile)
                sys.stderr.flush()

        return myspps
    
    def dataframe(self, boldout):
        
        df = {}
        for line in open(boldout, 'r').readlines():

            line = line.strip()
            seq,spps,perc,_ = line.split("\t")
            tmp_perc = float(perc)
            # {
            #  seq : {spps:val, spps2:val2},
            #  seq2: {spps:val, spps2:val2},
            # }
            if df.__contains__(seq): 

                if df[seq].__contains__(spps):
                    prior_per = df[seq][spps]

                    if tmp_perc > prior_per:
                        df[seq][spps] = tmp_perc
                else:
                    df[seq].update( { spps:tmp_perc} )
            else:
                df[seq] = { spps:tmp_perc }

        return df
    
    def classify(self, df, spps_aln):
        
        row    = []
        ROWSTR = "{seq}\t{species}\t{matchtype}\n"

        for k,v in spps_aln.items():

            if df.__contains__(k):
                spps_perc = dict( filter( lambda kv: kv[1] >= self.threshold, df[k].items() ) )

                if spps_perc:            
                    matched_spps = set(list(spps_perc))

                    if v in matched_spps:

                        if matched_spps.__len__() == 1:
                            row.append(
                                ROWSTR.format(
                                    seq = k,
                                    species =  matched_spps.pop(),
                                    matchtype = 'match'))
                        else:
                            row.append(
                                ROWSTR.format(
                                    seq = k,
                                    species = ','.join(matched_spps),
                                    matchtype = 'ambiguous'))
                    else:
                        row.append(
                            ROWSTR.format(
                                seq = k,
                                species = ','.join(matched_spps),
                                matchtype = 'mismatch'))         
                else:
                    bmatch = sorted(df[k], key = df[k].get, reverse = True)[0:3]
                    row.append(
                        ROWSTR.format(
                            seq = k,
                            species = ','.join(bmatch),
                            matchtype = 'below %s' % self.threshold))
            else:
                row.append(
                    ROWSTR.format(
                        seq = k,
                        species = 'NA',
                        matchtype = 'NA'))

        return row

    def id_engine(self):

        boldin  = self.sequence + "_bold_in"
        boldout = self.sequence + "_bold_out"
        classout = self.sequence + "_bold" if not self.outfile else self.outfile
        
        # species-level dictionary
        spps_aln  = self.checkspps()
        # whole fastas
        tmp_fasta = fas_to_dic(self.dealignment)        
        # trimmed fasta
        with open(boldin, 'w') as f:
            for i in spps_aln:
                f.write( ">%s\n%s\n" % ( i, tmp_fasta[">" + i] ) ) 
                
        id_engine(query       = boldin, 
                  db          = self.bolddatabase,
                  make_blast  = self.make_blast,
                  quiet       = self.quiet,
                  threshold   = self.threshold,
                  fileoutname = boldout)
        
        mydf       = self.dataframe(boldout)
        classified = self.classify(mydf, spps_aln)
        
        with open(classout, 'w') as f:
            f.write("sequence\tbest_match\tmatch_type\n")
            f.writelines(classified)
            
        if self.removeintermediate:            
            os.remove(self.dealignment)
            os.remove(boldin)
            os.remove(boldout)
            os.remove(boldout + "_filtered")
