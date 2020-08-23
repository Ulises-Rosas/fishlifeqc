import os
from multiprocessing import Pool
from fishlifeqc.utils import fas_to_dic

class Deletion:

    def __init__(self, 
                controlfile = None, 
                filetype    = 'rblast',
                threads     = 1):

        self.controlfile = controlfile
        self.filetype    = filetype
        self.suffix      = "_%sd" % filetype
        self.threads     = threads

        self.FAILEDFILE  = "failed_deletion.txt"

    def readcontrolfile(self):

        if self.filetype == 'rblast':

            files = {}

            for i in open(self.controlfile, 'r').readlines()[1:]:

                filen,sample,_ = i.split(',')

                if not files.__contains__(filen):
                    files[filen] = [sample]

                else:
                    files[filen] += [sample]

        return files

    def constraint(self, header_todel):

        myfile, headers = header_todel

        try:
            myfasta = fas_to_dic(myfile)

        except FileNotFoundError:
            return myfile

        myfasta  = fas_to_dic(myfile)

        newfasta = {}
        for h,s in myfasta.items():

            if not h.replace(">", "") in headers:
                newfasta[h] = s

        with open( myfile + self.suffix, 'w') as f:
            for h,s in newfasta.items():
                f.write("%s\n%s\n" % (h,s))

        return None

    def run(self):

        files = tuple(self.readcontrolfile().items())

        with Pool(processes = self.threads) as p:
            failed = [*p.map(self.constraint, files)]

        failed = list(filter(None, failed))

        if failed:
            with open(self.FAILEDFILE, 'w') as f:
                for i in failed:
                    f.write( i + "\n")
# Deletion(
#     controlfile= "mismatch_pairedblastn.txt",
#     filetype='rblast',
#     threads= 1
#     ).run()