import os
import sys
import csv
from multiprocessing import Pool
from fishlifeqc.utils import fas_to_dic

class Deletion:

    def __init__(self, 
                sequences   = None,
                controlfile = None,
                deletion_list = None,
                filetype    = 'list',
                threads     = 1):

        self.deletion_list = deletion_list # plain list of species
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

        if not self.sequences:
            return None

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

    def _read_deletion_list(self):
        out = []

        if not self.deletion_list:
            return out

        with open(self.deletion_list, 'r') as f:
            for i in f.readlines():
                line = i.strip()
                if line:
                    out.append(line)
        return out

    def delete_by_custom_list(self, aln: dict, custom_list: list = None):

        if not aln:
            return None

        if not custom_list:
            return aln

        out = {}
        for k,v in aln.items():
            header = k.strip().replace(">", "")
            if not header in custom_list:
                out[k] = v

        if not out:
            return None

        return out    
        
    def customCountTrimming(self, aln  : dict,
                            min_count  : int = 30, 
                            custom_list: list = None) -> dict:
        """
        Moves:
        1. delete by using a custom list (if not empty)
        2. delete by minimum amount of sequence headers

        returns None if not sequences
        """

        aln = self.delete_by_custom_list(aln, custom_list)

        if not aln:
            return None

        nheader = len(list(aln))

        if nheader < min_count:
            return None

        return aln

# self = Deletion(
#     sequences   =  ['../scripts/E0085.listd_allsets.NT_aligned_renamed.fasta_trimmed_round2_listdno_lava'],
#     controlfile =  "./../scripts/custom_deletion_list.txt",
#     filetype    = '_listd',
#     threads     =  1
# )
# """
# Microstomatidae_Nansenia_ardesiaca_EPLATE_46_F04
# Opisthoproctidae_Rhynchohyalus_natalensis_EPLATE_49_G08
# Galaxiidae_Galaxias_maculatus_EPLATE_25_B08
# Retropinnidae_Retropinna_tasmanica_EPLATE_47_A10
# Myxinidae_Myxine_knappi_EPLATE_24_D01
# """