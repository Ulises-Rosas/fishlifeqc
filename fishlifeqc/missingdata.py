
import csv
import sys
import fishlifeseq
import collections
from multiprocessing import Pool
from fishlifeqc.utils import isfasta, fas_to_dic, compare_alns
from fishlifeqc.delete import Deletion

FAILEDTOTRIM = "no_passed_trimming.txt"
MDATA_REPORT = "mdata_report.csv"

# source https://github.com/biopython/biopython/blob/master/Bio/Data/CodonTable.py
STOP_CODON_TABLE = {
    1:  {'codons': ['TAA', 'TAG', 'TGA'], 'name': 'Standard'},
    2:  {'codons': ['TAA', 'TAG', 'AGA', 'AGG'], 'name': 'Vertebrate Mitochondrial'},
    3:  {'codons': ['TAA', 'TAG'], 'name': 'Invertebrate Mitochondrial'},
    4:  {'codons': ['TAA', 'TAG', 'TGA'], 'name': 'Bacterial, Archaeal and Plant Plastid'},
    5:  {'codons': ['TAA', 'TAG'], 'name': 'Yeast Mitochondrial'},
    6:  {'codons': ['TAA', 'TAG'],'name': 'Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate '},
    7:  {'codons': ['TGA'], 'name': 'Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear'},
    8:  {'codons': ['TAA', 'TAG'], 'name': 'Echinoderm Mitochondrial; Flatworm Mitochondrial'},
    9:  {'codons': ['TAA', 'TAG'], 'name': 'Euplotid Nuclear'},
    10: {'codons': ['TAA', 'TAG', 'TGA'], 'name': 'Alternative Yeast Nuclear'},
    11: {'codons': ['TAA', 'TAG'], 'name': 'Ascidian Mitochondrial'},
    12: {'codons': ['TAG'], 'name': 'Alternative Flatworm Mitochondrial'},
    13: {'codons': ['TAA', 'TGA'], 'name': 'Blepharisma Macronuclear'},
    14: {'codons': ['TAA', 'TGA'], 'name': 'Chlorophycean Mitochondrial'},
    15: {'codons': ['TAA', 'TAG'], 'name': 'Trematode Mitochondrial'},
    16: {'codons': ['TCA', 'TAA', 'TGA'],'name': 'Scenedesmus obliquus Mitochondrial'},
    17: {'codons': ['TTA', 'TAA', 'TAG', 'TGA'], 'name': 'Thraustochytrium Mitochondrial'},
    18: {'codons': ['TAA', 'TAG'], 'name': 'Pterobranchia Mitochondrial'},
    19: {'codons': ['TAA', 'TAG'],'name': 'Candidate Division SR1 and Gracilibacteria'},
    20: {'codons': ['TAA', 'TAG', 'TGA'],'name': 'Pachysolen tannophilus Nuclear'},
    21: {'codons': ['TGA'], 'name': 'Karyorelict Nuclear'},
    22: {'codons': ['TAA', 'TAG', 'TGA'], 'name': 'Condylostoma Nuclear'},
    23: {'codons': ['TGA'], 'name': 'Mesodinium Nuclear'},
    24: {'codons': ['TGA'], 'name': 'Peritrich Nuclear'},
    25: {'codons': ['TAA', 'TAG'], 'name': 'Blastocrithidia Nuclear'},
    26: {'codons': ['TAA', 'TGA'], 'name': 'Balanophoraceae Plastid'},
    27: {'codons': ['TAG'], 'name': 'Cephalodiscidae Mitochondrial'}
}

class Missingdata(Deletion):

    def __init__(self,
                fastas = None, 
                htrim = 0.6,  # if not horizontal, just set self.htrim = 1
                vtrim = 0.1,
                itrim = 0.1,
                min_prop_sequences_per_aln = 0.25,
                min_prop_alns_per_sequence = 0,
                custom_deletion_list = None,
                outputsuffix = "_trimmed", 
                codon_aware = False,
                stop_opt    = 1,
                threads     = 1):

        # ------ share variables between objs -------

        self.deletion_list = custom_deletion_list # specific list of taxa to delete on each alignment
        self.customlist    = self._read_deletion_list()

        # ------ share variables between objs -------

        self.fasta = fastas
        self.htrim = htrim # maximum gap percentage allowed per sequence
        self.vtrim = vtrim # maximum gap percentage allowed per columm
        self.itrim = itrim # maximum gap percentage allowed per internal columm
        self.min_prop_sequences_per_aln = min_prop_sequences_per_aln # Minimum proportion of sequences per alignments
        self.min_prop_alns_per_sequence = min_prop_alns_per_sequence

        self.min_sequences_per_aln = 4 # safe internal default
        
        self.outputsuffix = outputsuffix


        self.codon_aware  = codon_aware
        #   |
        #    -> if codon_aware:
        self.stop_codon_table = STOP_CODON_TABLE
        self.stop_opt         = stop_opt

        self.threads = threads

        #  ------- DEPRECATED ------- #
        # stops for mitchondrial DNA
        self.mt_stop_codons=["TAA", "TAG", "AGA", "AGG"]
        # stops for standard genetic code
        self.std_stop_codons=["TAA", "TAG", "TGA"]
        #  ------- DEPRECATED ------- #

    @property
    def stop_codons(self):
        return self.stop_codon_table[self.stop_opt]['codons']

    def print_stop_lib(self):
        """
        print table of options
        for stop_lib
        """
        num_opts   = len( list(self.stop_codon_table) )
        name_span  = max([ len(v['name']) for v in self.stop_codon_table.values() ])
        codon_span = max([ len(", ".join(v['codons'])) for v in self.stop_codon_table.values() ])

        str_format = "%-7s| %-" + str(name_span) + "s| %-" + str(codon_span)+ "s"
        print()
        print( str_format   % ( "Option", "Name", "Codons" )           )
        print( "%s+-%s+-%s" % ("-"*7, "-"*name_span, "-"*codon_span )  )

        for opt in range(1, num_opts + 1):
            name      = self.stop_codon_table[opt]['name']
            codons    = self.stop_codon_table[opt]['codons']
            condons_f = ", ".join(codons)

            print( str_format % (opt, name, condons_f) )

        print()

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

    def horizontalTrimming(self, aln):
        """ 
        Moves:

        1. remove those sequences when
        the gap proportion is more
        than the threshold value

        2. Close only-gaps columns
        by above move

        * If `aln` is empty, returns None
        """
        if not aln:
            return None

        seqlength = Missingdata._seq_length(aln)

        out = {}
        for k,v in aln.items():
            if v.count('-')/seqlength <= self.htrim:
                out.update( {k:v}  )

        if not out:
            return None

        return Missingdata.close_gaps(out, self.codon_aware)

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
        return codon in self.stop_codons

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
        """
        - If rigth trimming reach left trimming, None is returned

        -  when more than two stop codons are present,
            `filter_stopcodons` returns an empty dict.
            Otherwise, it returns `aln` with lengths!

        - if `trimmed`, `self.close_gaps` closes 
            possible only-gaps columns 
            produced by `filter_stopcodons`
        """
        
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
            """   
            horizontal trimming (check stop codons).
            when more than two stop codons are present,
            `filter_stopcodons` returns an empty dict.
            Otherwise, it returns `aln` with lengths!
            """
            trimmed = self.filter_stopcodons(trimmed)
            """
            close possible only-gaps columns
            produced by `filter_stopcodons`
            """
            return Missingdata.close_gaps(aln = trimmed, is_codon_aware= True)
        else:
            return None

    @staticmethod
    def _seq_length(aln:dict)-> int:
        return len(next(iter(aln.values())))

    @staticmethod
    def _get_gap_prop(aln: dict, seqlength: int, offset: int)-> dict:

        idxs = []
        for c in range(0, seqlength, offset):

            mystr = ""
            for v in aln.values():
                mystr += v[c:c+offset]

            gap_prop = mystr.count('-')/len(mystr)
            idxs.append( (c, gap_prop) )

        return idxs

    @staticmethod
    def _filter(indexes: list, aln: dict, seqlength: int, offset: int) -> dict:
        """
        if indexes list is empty, it means
        there is anything to save because any column
        pass the maximun allowed proportion condition.

        then, return None
        """

        if not indexes:
            return None

        out = {}
        for k,v in aln.items():
            mystr  = ""
            for c in range(0, seqlength, offset):
                if c in indexes:
                    mystr += v[c:c+offset]
            out[k] = mystr

        seqlength = Missingdata._seq_length(out)

        if not seqlength:
            return None

        return out

    @staticmethod
    def internal_trimmer(aln: dict, is_codon_aware: bool = False, allowed_gap_prop: float = 0.1) -> dict:
        # aln = {'a':'---cat-aa', 'b':'---cataaa', 'c':'---cataaa'}
        """
        * If `aln` is empty, returns None

        * If not seqlength after processing, None is returned
        """

        if not aln:
            return None

        offset   = 3 if is_codon_aware else 1
        # seqlength = len( list( aln.values() )[0] )
        seqlength = Missingdata._seq_length(aln)

        idxs = []
        idxs_prop = Missingdata._get_gap_prop(aln, seqlength, offset)
        for c, gap_prop in idxs_prop:
            if gap_prop <= allowed_gap_prop:
                idxs.append(c)

        return Missingdata._filter(idxs, aln, seqlength, offset)

    @staticmethod
    def convertN2Gap(aln: dict) -> dict:
        """
        `Missingdata.convertN2Gap` is used as
        both aligner and `self.filter_stopcodons` add "N"s.
        Used once if called from other module (e.g.,`monophyly.py`)
        """

        for k,v in aln.items():
            aln[k] = v.replace("N", "-")

        return aln

    @staticmethod
    def close_gaps(aln: dict, is_codon_aware: bool = True) -> dict:
        """
        Moves:
        1. It converts N to gaps
        2. get proportions
        3. filter

        It closes the empty columns formed by:
        - `self.filter_stopcodons`
        - `self.horizontalTrimming`

        * If `aln` is empty, returns None
        * If `aln` is not empty, `aln` has values with length
        * If not seqlength after processing, None is returned
        """
        # aln = {'a':'---cat--a', 'b':'---cat--a'}
        # is_codon_aware = not True
        if not aln:
            return None
            
        offset    = 3 if is_codon_aware else 1
        aln       = Missingdata.convertN2Gap(aln = aln)
        seqlength = Missingdata._seq_length(aln)

        idxs = []
        idxs_prop = Missingdata._get_gap_prop(aln, seqlength, offset)

        for c, gap_prop in idxs_prop:
            if not gap_prop == 1:
                idxs.append(c)

        return Missingdata._filter(idxs, aln, seqlength, offset)
        
    def writeresults(self, obj, name):

        outname = name + self.outputsuffix

        with open(outname, 'w') as f:

            for k,v in obj.items():
                f.write( "%s\n%s\n" % (k,v))

    def verticalTrimming(self, aln: dict) -> dict:
        """
        Moves:
        1. Trim edges (None if all is removed)
        2. Trim internal (None if all is removed)

        returns either None or a filled dict
        """
        if not aln:
            return None

        if not self.codon_aware:
            # vertical trim, length change
            trimmed = self.trimedges(aln = aln)
        else:
            # vertical trim, length change &
            # horizontal trimming (check stop codons)
            trimmed = self.trimedges_codonaware(aln = aln)

        # None or filled dict
        trimmed = Missingdata.internal_trimmer(
                                    aln              = trimmed,
                                    is_codon_aware   = self.codon_aware,
                                    allowed_gap_prop = self.itrim)
        return trimmed
    
    def trimiterator(self, fasta):
        """
        Main trimmer
        ============

        Moves:

        1. Trimming by sequence count
            and custom list of taxa:
            minimum sequence count 
            per sequence is setted at
            25% of the number of all entered
            sequences by default and the 
            custom list is optional 
            (i.e., default is None)

        2. vertical trimming:
            - edges (codon aware)
            - internal colums

        3. horizontal trimming:
            Looking for individual 
            sequences with high proportion
            of gaps

        Sequences are pre-processed by
        converting N to "-" due to some aligners
        add "N"s to increase the score of the alignment.
        This is particularly true for the MACSE aligner
        whose algorithm can add up to two N when it cannot fit
        a base to any surrounding codons.
        """


        sys.stderr.write("Trimming: %s\n" % fasta )
        sys.stderr.flush()

        # fasta = self.fasta[2]
        if not isfasta(fasta):
            return None

        alignment = fas_to_dic(fasta)
        alignment = self.convertN2Gap(alignment)
        seqlength = self.check_aln(alignment, fasta)

        if not seqlength:
            return None

        # trim by count and a custom list
        trimmed = self.customCountTrimming(alignment, self.min_sequences_per_aln, self.customlist)
        # if not trimmed:
        #     return None

        # close possible only-gaps columns
        # produced by `customCountTrimming` 
        # (i.e., only-gap columns supported only by 
        #        sequences belonging to `self.customlist`)
        trimmed = self.close_gaps(aln = trimmed, is_codon_aware =  self.codon_aware)

        trimmed = self.verticalTrimming(aln = trimmed)

        # horizontal trim, not length change
        trimmed = self.horizontalTrimming(aln = trimmed)

        if trimmed:
            self.writeresults(trimmed, fasta)
            desc = compare_alns(alignment, trimmed)

            return (fasta, desc)
        else:
            return None

    def _update_custom_list(self, headers_f_joined):

        headers_count = collections.Counter(headers_f_joined)
        thresh = round( len(self.fasta) * self.min_prop_alns_per_sequence )
        to_del = [k.strip().replace(">", "") for k,v in headers_count.items() if v < thresh]

        self.customlist += to_del

    def _update_min_counts(self, all_headers):
        """
        set the minimum amount of sequences
        that each alignment must have
        this minimum values obtained after
        get the total number  from all alignments used as
        input
        """
        headers_f = list( filter(None, all_headers) )

        if not headers_f:
            return None

        headers_f_joined = []
        for h in headers_f:
            headers_f_joined.extend(h)
        
        headers_uniq = set(headers_f_joined)
        new_min_seqs_per_aln = round( len(headers_uniq) * self.min_prop_sequences_per_aln )

        if new_min_seqs_per_aln > self.min_sequences_per_aln:
            self.min_sequences_per_aln = new_min_seqs_per_aln

        self._update_custom_list(headers_f_joined)

    def _write_report(self, table, failed):
        # table = descriptions
        table_f = list(filter(None, table))

        if not table_f:
            sys.stderr.write("\nAny alignment passed trimming steps\n")
            sys.stderr.flush()
            return None

        exons_before = len(self.fasta)
        exons_after  = exons_before - len(failed)
   
        site_before = 0
        gaps_before = 0
        headers_before = 0

        site_after = 0
        gaps_after = 0
        headers_after = 0 

        for line in table_f:
            site_before += line[0]
            gaps_before += line[1]
            headers_before += line[2]

            site_after += line[3]
            gaps_after += line[4]
            headers_after += line[5]

        exon_perc    = round(exons_after*100/exons_before, 2)
        sites_perc   = round(site_after*100/site_before, 2)
        gaps_perc    = round(gaps_after*100/gaps_before, 2)
        headers_perc = round(headers_after*100/headers_before, 2)

        out  = [[ 'description', 'before',       'after',        'percentage'],
                [ 'alignments' , exons_before,   exons_after,    exon_perc   ],
                [ 'sites'      , site_before,    site_after,     sites_perc  ],
                [ 'gaps'       , gaps_before,    gaps_after,     gaps_perc   ],
                [ 'headers'    , headers_before, headers_after,  headers_perc]]

        with open( MDATA_REPORT, 'w') as f:
            writer = csv.writer(f)
            writer.writerows( out )

        sys.stderr.write("\nReport written at %s\n" % MDATA_REPORT)
        sys.stderr.flush()

    def run(self):

        with Pool(processes = self.threads) as p:

            if self.min_prop_sequences_per_aln or self.min_prop_alns_per_sequence:
                sys.stderr.write("\nMapping taxa...\r")
                sys.stderr.flush()

                all_headers = [*p.map(fishlifeseq.headers, self.fasta)]
                self._update_min_counts(all_headers)

                sys.stderr.write("Mapping taxa...Ok\n\n")
                sys.stderr.flush()

            descriptions = []
            passed       = []

            # for fasta in self.fasta:
            #     preout = p.map_async(self.trimiterator, (fasta,)).get()[0]
            #     if preout:
            #         fasta, desc = preout
            #         passed.append(fasta)
            #         descriptions.append(desc)

            preout = []
            for fasta in self.fasta:
                result  = p.map_async(self.trimiterator, (fasta,))
                preout.append(result)

            for pr in preout:
                gotit = pr.get()[0]
                if gotit:
                    fasta, desc = gotit
                    passed.append(fasta)
                    descriptions.append(desc)

        failed = set(self.fasta) - set(list(filter(None, passed)))

        if failed:
            with open(FAILEDTOTRIM, 'w') as f:
                for i in failed:
                    f.write( "%s\n" % i)

        self._write_report(descriptions, failed)

# --- testings -- #
# import glob
# fastas = glob.glob("/Users/ulises/Desktop/test/mdata/*.fasta")
# fastas = [
#     '/Users/ulises/Desktop/test/mdata/E1532.fasta',
#     '/Users/ulises/Desktop/test/mdata/E0219.fasta',
#     '/Users/ulises/Desktop/test/mdata/E1078.fasta'
#  ]
# # # fishlifeseq.headers(fastas[0])
# custom_deletion_list = "/Users/ulises/Desktop/test/mdata/custom_deletion_list.txt"
# self = Missingdata(
#     custom_deletion_list=custom_deletion_list,
#     fastas=fastas,
#     outputsuffix = "_trimmed.fasta",
#     itrim= 0.6,
#     vtrim= 0.5,
#     min_prop_sequences_per_aln = 0.25,
#     min_prop_alns_per_sequence = 0.25,
#     codon_aware=True,
#     threads= 3
# )
# self.print_stop_lib()
# print(self.min_count)
# -- testings -- #