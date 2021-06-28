
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
                isaminoacid = False,
                write_bad   = False,
                trim_seqs   = False,
                pval        = 0.05,
                threads     = 1):

        self.sequences   = sequences
        self.symtype     = symtype
        self.nexusformat = nexusformat
        self.isaminoacid = isaminoacid
        self.write_bad   = write_bad
        self.trim_seqs   = trim_seqs
        self.pval        = pval
        self.codon_aware = codon_aware

        self.suffix    = suffix
        self.notpassed = 'NotPassed_' + suffix + ".txt"
        self.threads   = threads


        self.trim_nex_format = "\tcharset {pos} = {file}: {pos_num}-{seq_len}\\{npart};\n"
        self.trim_rax_format = "DNA, {file}_{pos} = {pos_num}-{seq_len}\\{npart}\n"

        # knock down variables
        self.replace_char = 'N'
        self.controlfile = None
        self.possible_codons = ['pos1', 
                                'pos2', 
                                'pos3']



    def _write_fasta(self, aln, aln_file, pval):
        # aln_base = os.path.basename(aln_file)
        if pval and pval >= self.pval:

            outfile = "%s_%s" % (aln_file, self.suffix)
            with open(outfile, 'w') as f:
                for k,v in aln.items():
                    f.write( "%s\n%s\n" % (k,v) )

            return True

        else:
            if self.write_bad:
                outfile_np = "%s_%s_NotPassed" % (aln_file, self.suffix)
                # outfile = "%s_%s" % (aln_file, self.suffix)
                with open(outfile_np, 'w') as f:
                    for k,v in aln.items():
                        f.write( "%s\n%s\n" % (k,v) )

            return False

    def trim_columns(self, to_rm, aln, seq_len):

        if to_rm:
            trimmed = {}
            for k,v in aln.items():
                mystr = ""
                for i in range(0, seq_len, 3):

                    F = '' if 'pos1' in to_rm else v[i]
                    S = '' if 'pos2' in to_rm else v[i + 1]
                    T = '' if 'pos3' in to_rm else v[i + 2]
                    mystr += (F + S + T)

                trimmed[k] = mystr

            trimmed_len = len(next(iter(trimmed.values())))

            return (trimmed, trimmed_len)

        else:
            return (aln, seq_len)

    def set_format(self, pos_num, file, seq_len, npart):

        if self.nexusformat:

            return self.trim_nex_format.format(
                        pos     = "pos%s" % pos_num,
                        file    = os.path.basename(file),
                        pos_num = pos_num,
                        seq_len = seq_len,
                        npart   = npart
                    )

        else:
            return self.trim_rax_format.format(
                        pos     = "pos%s" % pos_num,
                        file    = os.path.basename(file),
                        pos_num = pos_num,
                        seq_len = seq_len,
                        npart   = npart
                    )

    def _create_part_and_trim(self, another_pval_table, seq_len, aln, aln_file):

        passed_pos = []
        not_passed_pos = []

        outfile = "%s_%s" % (aln_file, self.suffix)
        outfile_np = "%s_%s_NotPassed" % (aln_file, self.suffix)

        for pos, pval in another_pval_table.items():

            if pval and pval >= self.pval:
                passed_pos.append(pos)

            else:
                not_passed_pos.append(pos)

        passed_aln, passed_aln_len = self.trim_columns(not_passed_pos, aln, seq_len)
        not_passed_aln, not_passed_aln_len = self.trim_columns(passed_pos, aln, seq_len)

        if passed_pos:
            part_file = []
            npart = len(passed_pos)

            for i in range(npart):

                part_file.append( 
                    self.set_format(
                        pos_num = i + 1,
                        file    = outfile, 
                        seq_len = passed_aln_len, 
                        npart   = npart
                    )
                )

            if self.nexusformat:
                part_file = ["#nexus\n", "begin sets;\n"] + part_file + ['end;\n']

            with open( outfile + ".nex", 'w' ) as f:
                f.writelines(part_file)

            with open( outfile, 'w' ) as f:
                for k,v in passed_aln.items():
                    f.write( "%s\n%s\n" % (k,v) )

        if self.write_bad and not_passed_pos:

            part_file = []
            npart = len(not_passed_pos)

            for i in range(npart):

                part_file.append( 
                    self.set_format(
                        pos_num = i + 1,
                        file    = outfile_np, 
                        seq_len = not_passed_aln_len, 
                        npart   = npart
                    )
                )

            if self.nexusformat:
                part_file = ["#nexus\n", "begin sets;\n"] + part_file + ['end;\n']

            with open( outfile_np + ".nex", 'w' ) as f:
                f.writelines(part_file)

            with open( outfile_np, 'w' ) as f:
                for k,v in not_passed_aln.items():
                    f.write( "%s\n%s\n" % (k,v) )

        return not_passed_pos


    def _create_part_files(self, another_pval_table, seq_len, aln, aln_file):

        outfile    = "%s_%s" % (aln_file, self.suffix)
        outfile_np = "%s_%s_NotPassed" % (aln_file, self.suffix)

        # outfile_base = 'ads'
        outfile_base = os.path.basename(outfile)
        outfile_np_base = os.path.basename(outfile_np)

        nex_format = "\tcharset {pos} = {file}: {pos_num}-{seq_len}\\3;\n"
        rax_format = "DNA, {file}_{pos} = {pos_num}-{seq_len}\\3\n"

        _passed_part = []
        passed_pos = []
        _not_passed_part = []
        not_passed_pos = []

        for pos, pval in another_pval_table.items():

            pos_num = pos.replace("pos", "")
            
            if pval and pval >= self.pval:

                if self.nexusformat:
                    _passed_part.append( 
                        nex_format.format(
                                pos = pos, 
                                pos_num = pos_num,
                                file = outfile_base, 
                                seq_len = seq_len
                        )
                    )
                    
                else:
                    _passed_part.append(
                        rax_format.format(
                                pos = pos, 
                                pos_num = pos_num,
                                file = outfile_base, 
                                seq_len = seq_len
                            )
                    )
                passed_pos.append(pos)

            else:
                if self.nexusformat:
                    _not_passed_part.append(
                        nex_format.format(
                            pos = pos, 
                            pos_num = pos_num,
                            file = outfile_np_base, 
                            seq_len = seq_len
                            )
                        )

                else:
                    _not_passed_part.append(
                        rax_format.format(
                                pos = pos, 
                                pos_num = pos_num,
                                file = outfile_np_base, 
                                seq_len = seq_len
                                )
                        )

                not_passed_pos.append(pos)

        if self.nexusformat:
            if _passed_part:
                _passed_part = ["#nexus\n", "begin sets;\n"] + _passed_part + ['end;\n']

            if _not_passed_part:
                _not_passed_part = ["#nexus\n", "begin sets;\n"] + _not_passed_part + ['end;\n']

        if _passed_part:
            with open( outfile, 'w' ) as f:
                for k,v in aln.items():
                    f.write( "%s\n%s\n" % (k,v) )

            with open( outfile + ".nex", 'w' ) as f:
                f.writelines(_passed_part)

        if self.write_bad and _not_passed_part:

            with open( outfile_np, 'w' ) as f:
                for k,v in aln.items():
                    f.write( "%s\n%s\n" % (k,v) )

            with open( outfile_np + ".nex", 'w' ) as f:
                f.writelines(_not_passed_part)

        return not_passed_pos
  
    def _write_fasta_c(self, aln, seq_len, aln_file, symtype, pval_table):

        another_pval_table = {
            'pos1': pval_table['pos1'][symtype],
            'pos2': pval_table['pos2'][symtype],
            'pos3': pval_table['pos3'][symtype],
        }

        if self.trim_seqs:
            return self._create_part_and_trim( 
                        another_pval_table = another_pval_table,
                        seq_len  = seq_len, 
                        aln      = aln,
                        aln_file = aln_file
                    )

        else:
            return self._create_part_files(             
                            another_pval_table = another_pval_table, 
                            seq_len  = seq_len,
                            aln      = aln,
                            aln_file = aln_file,
                        )

    def sym_iterator(self, aln_file):
        # aln_file = "/Users/ulises/Desktop/GOL/data/alldatasets/nt_aln/internally_trimmed/malns_36_mseqs_27/round2/no_lavaretus/protein/round1/E0109.aaseqs.fasta"
        
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
             IntPval) = self._symmetries(all_pairs, aln, aa = self.isaminoacid)

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

            if not (seq_len % 3):
                sys.stderr.write("'%s' sequences length not divisible by three" % aln_base)
                sys.stderr.flush()
                return None

            codon1, codon2, codon3 = self._split_aln(aln, seq_len)
            
            (SymPval_pos1,
             MarPval_pos1,
             IntPval_pos1) = self._symmetries(all_pairs, codon1, aa = self.isaminoacid)

            (SymPval_pos2,
             MarPval_pos2,
             IntPval_pos2) = self._symmetries(all_pairs, codon2, aa = self.isaminoacid)

            (SymPval_pos3,
             MarPval_pos3,
             IntPval_pos3) = self._symmetries(all_pairs, codon3, aa = self.isaminoacid)

            pval_table = {
                'pos1': {
                    'sym': SymPval_pos1,
                    'mar': MarPval_pos1,
                    'int': IntPval_pos1,
                },

                'pos2': {
                    'sym': SymPval_pos2,
                    'mar': MarPval_pos2,
                    'int': IntPval_pos2,
                },

                'pos3': {
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

    def knockdown_iterator(self, aln_file):
        # aln_file = "/Users/ulises/Desktop/GOL/software/fishlifeqc/demo/para/E1381.fasta"
        aln_base = os.path.basename(aln_file)
        outfile  = "%s_%s" % (aln_file, self.suffix)

        sys.stdout.write("Processing: %s\n" % aln_base)
        sys.stdout.flush()

        aln = fas_to_dic(aln_file)

        if not self._control_file.__contains__(aln_file):

            with open(outfile, 'w') as f:
                for k,v in aln.items():
                    f.write("%s\n%s\n" % (k,v)) 

            return None

        pos_rm = self._control_file[aln_base]

        # aln_file = os.path.join(self.path, aln_base)
        # aln = fas_to_dic(aln_file)
        # os.rename(aln_file, "%s_%s" % (aln_file, self.suffix) )

        if len(pos_rm) == 3:
            return None

        seq_len = len(next(iter(aln.values())))
        trimmed = {}
        for k,v in aln.items():

            mystr = ""
            for i in range(0, seq_len, 3):

                F = self.r_char if 'pos1' in pos_rm else v[i]
                S = self.r_char if 'pos2' in pos_rm else v[i + 1]
                T = self.r_char if 'pos3' in pos_rm else v[i + 2]
                mystr += (F + S + T)

            trimmed[k] = mystr

        with open(outfile, 'w') as f:
            for k,v in trimmed.items():
                f.write("%s\n%s\n" % (k,v))

    def knockdown_columns(self, controlfile = None, replace_char = 'N'):

        self.r_char = replace_char
        self.controlfile = controlfile
        # self.controlfile =  "/Users/ulises/Desktop/GOL/software/\
        # fishlifeqc/demo/para/NotPassed_SymTest.txt"

        with Pool(processes = self.threads) as p:

            preout = []
            for fa in self.sequences:
                result  = p.map_async(self.knockdown_iterator, (fa,))
                preout.append(result)

            for pr in preout:
                pr.get()

# self = SymTests(symtype='sym', write_bad=True)


# not_passed = '/Users/ulises/Desktop/GOL/data/alldatasets/nt_aln/internally_trimmed/malns_36_mseqs_27/round2/no_lavaretus/srh_tests/NotPassed_SymTest.txt'
# partitions = '/Users/ulises/Desktop/GOL/data/alldatasets/nt_aln/internally_trimmed/malns_36_mseqs_27/round2/no_lavaretus/srh_tests/mypartitions.txt'

# not_passed_rows = []
# with open(not_passed, 'r') as f:
#     reader = csv.reader(f, delimiter = ",")
#     for row in reader:
#         not_passed_rows.append(row)

# not_passed_formated = ["%s_%s" % (i[0], i[1])  for i in not_passed_rows ]

# import re

# good_partitions = []
# bad_partitions = []
# with open(partitions, 'r') as f:
#     for line in f.readlines():

#         if not re.findall("=", line):
#             good_partitions.append(line)
#             bad_partitions.append(line)
#         else:
#             tmp_part = re.split("[ ]+", line)[2]

#             if tmp_part in not_passed_formated:
#                 bad_partitions.append(line)
#             else:
#                 good_partitions.append(line)


# with open(partitions + "_bad_ones", 'w') as f:
#     f.writelines(bad_partitions)

# with open(partitions + "_good_ones", 'w') as f:
#     f.writelines(good_partitions)
