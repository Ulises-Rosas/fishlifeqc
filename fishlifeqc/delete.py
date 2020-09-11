import os
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
        self.suffix      = "_%sd" % filetype
        self.threads     = threads

        # self.FAILEDFILE  = "failed_deletion.txt"
        self.customlist = []

    def readcontrolfile(self):

        if self.filetype == 'list':

            self.customlist = [i.strip() for i in open(self.controlfile, 'r').readlines()]

    def prot_headers(self, file):
        # file = './../data/mock1.txt'
        aln  = fas_to_dic(file)

        with open( file + self.suffix, 'w') as f:

            for k,v in aln.items():

                if not k.replace(">", "") in self.customlist:
                    f.write("%s\n%s\n" % (k,v))

    def headers(self):

        self.readcontrolfile()

        with Pool(processes = self.threads) as p:

            [*p.map(self.prot_headers, self.sequences)]


# Deletion(
#     sequences   =  ['data/COI.NT_aligned.fasta_trimmed',
#                     'data/E0537.NT_aligned.fasta_trimmed',
#                     'data/E1718.NT_aligned.fasta_trimmed'],
#     controlfile =   "data/mydellist.txt",
#     filetype    = 'list',
#     threads     =  2
# ).headers()
# """
# Microstomatidae_Nansenia_ardesiaca_EPLATE_46_F04
# Opisthoproctidae_Rhynchohyalus_natalensis_EPLATE_49_G08
# Galaxiidae_Galaxias_maculatus_EPLATE_25_B08
# Retropinnidae_Retropinna_tasmanica_EPLATE_47_A10
# Myxinidae_Myxine_knappi_EPLATE_24_D01
# """