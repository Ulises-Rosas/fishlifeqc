import os
import sys
import csv
from multiprocessing import Pool
from fishlifeqc.utils import fas_to_dic

class Deletion:

    def __init__(self, 
                sequences   = None,
                controlfile = None, 
                filetype    = 'list',
                threads     = 1):

        self.sequences   = sequences
        self.controlfile = controlfile
        self.filetype    = filetype
        self.suffix      = filetype
        self.threads     = threads

        # self.FAILEDFILE  = "failed_deletion.txt"
        self.customlist = []

    @property
    def _control_file(self):
        df = {}
        myrows = []

        with open(self.controlfile, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                myrows.append(row)

        for i in myrows:
            exon = i[0]
            header = i[1]
            if not df.__contains__(exon):
                df[exon] = [header]
            else:
                df[exon] += [header]

        return df

    def readcontrolfile(self):

        if self.filetype == 'list':

            self.customlist = [i.strip() for i in open(self.controlfile, 'r').readlines()]

    def prot_headers(self, file):

        if not os.path.exists(file):
            sys.stderr.write("\nFile '%s' does not exist\n" % file)
            sys.stderr.flush()
            return None

        aln  = fas_to_dic(file)

        with open( file + self.suffix, 'w') as f:

            for k,v in aln.items():

                if not k.replace(">", "") in self.customlist:
                    f.write("%s\n%s\n" % (k,v))

    def headers(self):

        self.readcontrolfile()

        with Pool(processes = self.threads) as p:

            [*p.map(self.prot_headers, self.sequences)]

    def _header_exon(self, seq):

        self.customlist = self._control_file[seq]
        self.prot_headers(seq)
        
    def header_exon(self):

        with Pool(processes = self.threads) as p:
            for seq in list(self._control_file):
                p.apply_async( self._header_exon, (seq,) ).get()

# Deletion(
#     sequences   =  ['data/COI.NT_aligned.fasta_trimmed',
#                     'data/E0537.NT_aligned.fasta_trimmed',
#                     'data/E1718.NT_aligned.fasta_trimmed'],
#     controlfile =  "head_headers_per_exon.csv",
#     filetype    = 'list',
#     threads     =  2
# ).header_exon()
# """
# Microstomatidae_Nansenia_ardesiaca_EPLATE_46_F04
# Opisthoproctidae_Rhynchohyalus_natalensis_EPLATE_49_G08
# Galaxiidae_Galaxias_maculatus_EPLATE_25_B08
# Retropinnidae_Retropinna_tasmanica_EPLATE_47_A10
# Myxinidae_Myxine_knappi_EPLATE_24_D01
# """