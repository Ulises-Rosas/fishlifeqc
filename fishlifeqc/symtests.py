
import os
import sys
import csv
import itertools
from multiprocessing import Pool

from fishlifeqc.utils import fas_to_dic, remove_files
from fishlifetraits.stats import Features



class SymTests(Features):
    def __init__(self, 
                sequences   = None,
                codon_aware = False,
                suffix      = 'SymTest',
                symtype     = 'sym',
                nexusformat = True,
                pval        = 0.05,
                threads     = 1):

        self.sequences   = sequences
        self.symtype     = symtype
        self.nexusformat = nexusformat
        self.pval        = pval
        self.codon_aware = codon_aware

        self.suffix    = suffix
        self.notpassed = 'NotPassed_' + suffix + ".txt"
        self.threads   = threads


        # knock down variables
        self.path = '.'
        self.controlfile = None
        self.possible_codons = ['pos_1', 
                                'pos_2', 
                                'pos_3']



    def _write_fasta(self, aln, aln_file, pval):
        # aln_base = os.path.basename(aln_file)
        if pval and pval >= self.pval:

            outfile = "%s_%s" % (aln_file, self.suffix)
            with open(outfile, 'w') as f:
                for k,v in aln.items():
                    f.write( "%s\n%s\n" % (k,v) )

            return True

        else:
            return False

    def _write_fasta_c(self, aln, seq_len, aln_file, symtype, pval_table):

        aln_base = os.path.basename(aln_file)
        outfile  = "%s_%s" % (aln_file, self.suffix)

        pval_c1 = pval_table['pos_1'][symtype] 
        pval_c2 = pval_table['pos_2'][symtype] 
        pval_c3 = pval_table['pos_3'][symtype] 

        part_file  = []
        not_passed = []

        if pval_c1 and pval_c1 >= self.pval:

            if self.nexusformat:
                part_file.append("\tcharset pos_1 = %s: 1-%s\\3;\n" % (aln_base, seq_len))

            else:
                part_file.append("DNA, %s_pos1 = 1-%s\\3\n" % (aln_base, seq_len))        
        else:
            not_passed.append('pos_1')

        if pval_c2 and pval_c2 >= self.pval:

            if self.nexusformat:
                part_file.append("\tcharset pos_2 = %s: 2-%s\\3;\n" % (aln_base, seq_len))

            else:
                part_file.append("DNA, %s_pos2 = 2-%s\\3\n" % (aln_base, seq_len))
        else:
            not_passed.append('pos_2')

        if pval_c3 and pval_c3 >= self.pval:

            if self.nexusformat:
                part_file.append("\tcharset pos_3 = %s: 3-%s\\3;\n" % (aln_base, seq_len))

            else:
                part_file.append("DNA, %s_pos3 = 3-%s\\3\n" % (aln_base, seq_len))
        else:
            not_passed.append('pos_3')

        if part_file:

            if self.nexusformat:
                part_file = ["#nexus\n", "begin sets;\n"] + part_file + ['end;\n']
                
            with open( outfile, 'w' ) as f:
                for k,v in aln.items():
                    f.write( "%s\n%s\n" % (k,v) )

            with open( outfile + ".nex", 'w' ) as f:
                f.writelines(part_file)

        return not_passed

    def sym_iterator(self, aln_file):
        # aln_file = "/Users/ulises/Desktop/GOL/data/alldatasets/nt_aln\
        # /internally_trimmed/malns_36_mseqs_27/\
        # round2/no_lavaretus/McL/E0795.r2_para_no_lavaretus_TBL_tlike_aln"
        
        aln_base = os.path.basename(aln_file)
        sys.stdout.write("Processing: %s\n" % aln_base)
        sys.stdout.flush()

        aln = fas_to_dic(aln_file)

        all_pairs = []
        for h1,h2 in itertools.combinations(aln, 2):
            all_pairs.append(  (h1.replace(">", ""), h2.replace(">", "")) )

        if not self.codon_aware:
            (SymPval,
             MarPval,
             IntPval) = self._symmetries(all_pairs, aln)

            if self.symtype == 'sym':
                passed = self._write_fasta(aln, aln_file, SymPval)

            elif self.symtype == 'mar':
                passed = self._write_fasta(aln, aln_file, MarPval)

            elif self.symtype == 'int':
                passed = self._write_fasta(aln, aln_file, IntPval)

            if not passed:
                with open( self.notpassed, 'a' ) as f:
                    writer = csv.writer(f, delimiter = ",")
                    writer.writerows( [[ aln_base, SymPval, MarPval, IntPval ]] ) 

        else:
            seq_len = len( next(iter(aln.values())) )
            codon1, codon2, codon3 = self._split_aln(aln, seq_len)
            
            (SymPval_pos1,
             MarPval_pos1,
             IntPval_pos1) = self._symmetries(all_pairs, codon1)

            (SymPval_pos2,
             MarPval_pos2,
             IntPval_pos2) = self._symmetries(all_pairs, codon2)

            (SymPval_pos3,
             MarPval_pos3,
             IntPval_pos3) = self._symmetries(all_pairs, codon3)

            pval_table = {
                'pos_1': {
                    'sym': SymPval_pos1,
                    'mar': MarPval_pos1,
                    'int': IntPval_pos1,
                },

                'pos_2': {
                    'sym': SymPval_pos2,
                    'mar': MarPval_pos2,
                    'int': IntPval_pos2,
                },

                'pos_3': {
                    'sym': SymPval_pos3,
                    'mar': MarPval_pos3,
                    'int': IntPval_pos3,
                }
            }

            not_passed = self._write_fasta_c(
                                aln      = aln,
                                seq_len  = seq_len, 
                                aln_file = aln_file, 
                                symtype  = self.symtype, 
                                pval_table = pval_table
                            )
 
            if not_passed:

                not_passed_table = []
                for pos in not_passed:

                    not_passed_table.append([
                        aln_base, 
                        pos, 
                        pval_table[pos]['sym'], 
                        pval_table[pos]['mar'], 
                        pval_table[pos]['int']
                    ])

                with open( self.notpassed, 'a' ) as f:
                    writer = csv.writer(f, delimiter = ",")
                    writer.writerows( not_passed_table ) 

    def init_files(self):

        if self.codon_aware:
            col_names = [["alignment", "position", "SymPval", "MarPval", "IntPval"]]
        else:
            col_names = [["alignment", "SymPval", "MarPval", "IntPval"]]

        with open(self.notpassed, "w") as f:
            writer = csv.writer(f, delimiter = ",")
            writer.writerows(col_names)

    def close_files(self):
        files_check = [
            self.notpassed
        ]

        to_rm = []

        sys.stdout.write("\n\n")
        sys.stdout.flush()

        for to_c in files_check:
            num_lines = sum( 1 for _ in open(to_c) )
            if num_lines == 1:
                to_rm.append(to_c)
            else:
                sys.stdout.write("'%s' report file was written\n" % to_c)
                sys.stdout.flush()
                
        remove_files(to_rm)

    def main(self):

        self.init_files()

        with Pool(processes = self.threads) as p:

            preout = []
            for fa in self.sequences:
                result  = p.map_async(self.sym_iterator, (fa,))
                preout.append(result)

            for pr in preout:
                pr.get()

        self.close_files()

    @property
    def _control_file(self):

        if not self.controlfile:
            return None

        # self.controlfile =  "/Users/ulises/Desktop/GOL/software/\
        # fishlifeqc/demo/para/NotPassed_SymTest.txt"
        df = {}
        myrows = []

        with open(self.controlfile, 'r') as f:
            reader = csv.reader(f, delimiter = ",")
            for row in reader:
                myrows.append(row)

        for i in myrows:
            exon = i[0]
            position = i[1]

            if not position in self.possible_codons:
                continue

            if not df.__contains__(exon):
                df[exon] = [position]

            else:
                df[exon] += [position]

        return df    

    def knockdown_iterator(self, aln_base):
        # aln_file = "/Users/ulises/Desktop/GOL/software/fishlifeqc/demo/para/E1381.fasta"
        sys.stdout.write("Processing: %s\n" % aln_base)
        sys.stdout.flush()

        pos_rm = self._control_file[aln_base]

        aln_file = os.path.join(self.path, aln_base)
        aln = fas_to_dic(aln_file)
        seq_len = len(next(iter(aln.values())))

        os.rename(aln_file, "%s_%s" % (aln_file, self.suffix) )

        if len(pos_rm) == 3:
            return None

        out = {}
        for k,v in aln.items():

            mystr = ""
            for i in range(0, seq_len, 3):

                F = '-' if 'pos_1' in pos_rm else v[i]
                S = '-' if 'pos_2' in pos_rm else v[i + 1]
                T = '-' if 'pos_3' in pos_rm else v[i + 2]
                mystr += (F + S + T)

            out[k] = mystr

        with open(aln_file, 'w') as f:
            for k,v in out.items():
                f.write("%s\n%s\n" % (k,v))

    def knockdown_columns(self, controlfile = None, path = "."):

        self.path = path
        self.controlfile = controlfile
        # self.controlfile =  "/Users/ulises/Desktop/GOL/software/\
        # fishlifeqc/demo/para/NotPassed_SymTest.txt"

        with Pool(processes = self.threads) as p:

            preout = []
            for fa in list(self._control_file):
                result  = p.map_async(self.knockdown_iterator, (fa,))
                preout.append(result)

            for pr in preout:
                pr.get()

# self = SymTests(symtype='sym')


