import os
import re
from multiprocessing import Pool
from fishlifeqc.utils import fas_to_dic

class Codonm:
    def __init__(self,
                 threads = 1,
                 alns    = None,
                 suffix  = ""):

        self.threads     = threads
        self.alns        = alns
        self.vercompname = "vertical_composition.txt" + suffix
        self.horcompname = "horizontal_composition.txt" + suffix

    def horizontal_composition(self, file):
        # file = "../data/E1718.NT_aligned.fasta_trimmed"

        basefile   = os.path.basename(file)
        aln        = fas_to_dic(file)

        seqlength  = len(aln[list(aln)[0]])
        lengthpercodon = seqlength/3

        nochar  = '(-|N)'
        binchar = re.compile(nochar)

        out     = []
        OUTROW  = "{exon},{sequence},{base},{pos1},{pos2},{pos3}\n"
        for header,sequence in aln.items():
            # header = '>Galaxiidae_Galaxias_truttaceus_EPLATE_25_B09'
            # sequence = aln[header]

            counter = {
                'pos1' : {'A': 0,
                          'G': 0,
                          'C': 0,
                          'T': 0,
                          'X': 0},
                'pos2' : {'A': 0,
                          'G': 0,
                          'C': 0,
                          'T': 0,
                          'X': 0},
                'pos3' : {'A': 0,
                          'G': 0,
                          'C': 0,
                          'T': 0,
                          'X': 0}
                }

            for i in range(0, seqlength, 3):

                F = sequence[i].upper()
                S = sequence[i + 1].upper()
                T = sequence[i + 2].upper()

                if binchar.match(F):
                    counter['pos1']['X'] += 1
                else:
                    counter['pos1'][F] += 1

                if binchar.match(S):
                    counter['pos2']['X'] += 1
                else:
                    counter['pos2'][S] += 1

                if binchar.match(T):
                    counter['pos3']['X'] += 1
                else:
                    counter['pos3'][T] += 1

            Om1 = lengthpercodon - counter['pos1']['X']
            Om2 = lengthpercodon - counter['pos2']['X']
            Om3 = lengthpercodon - counter['pos3']['X']
            
            for base in ['A', 'C', 'G', 'T']:

                out.append(
                    OUTROW.format(
                        exon     = basefile,
                        sequence = header.replace(">", ""),
                        base     = base,
                        pos1     = round(counter['pos1'][base]/Om1, 6) if Om1 else 0,
                        pos2     = round(counter['pos2'][base]/Om2, 6) if Om2 else 0,
                        pos3     = round(counter['pos3'][base]/Om3, 6) if Om3 else 0
                    )
                )

        return out

    def base_probability(self, column):
        # column  = F
        filterC = [i for i in column if re.findall("(A|C|G|T)", i)]

        # sample space
        Om = len(filterC)

        A = round(filterC.count('A')/Om, 6) if Om else 0
        C = round(filterC.count('C')/Om, 6) if Om else 0
        G = round(filterC.count('G')/Om, 6) if Om else 0
        T = round(filterC.count('T')/Om, 6) if Om else 0

        return (A,C,G,T)

    def vertical_composition(self, file):
        # file = "../data/E1718.NT_aligned.fasta_trimmed"
        basefile = os.path.basename(file)

        aln = fas_to_dic(file)

        seqlength    = len(aln[list(aln)[0]])
        counter = {
            'pos1' : {'A': [],
                    'G': [],
                    'C': [],
                    'T': []},
            'pos2' : {'A': [],
                    'G': [],
                    'C': [],
                    'T': []},
            'pos3' : {'A': [],
                    'G': [],
                    'C': [],
                    'T': []}
            }

        # iteration by codon
        for i in range(0, seqlength, 3):

            F = [v[i].upper()     for _,v in aln.items()]
            S = [v[i + 1].upper() for _,v in aln.items()]
            T = [v[i + 2].upper() for _,v in aln.items()]

            for pos,Co in [('pos1',F), ('pos2',S),('pos3',T)]:

                tmp_a,tmp_c,tmp_g,tmp_t = self.base_probability(Co)

                counter[pos]['A'] += [tmp_a]
                counter[pos]['C'] += [tmp_c]
                counter[pos]['G'] += [tmp_g]
                counter[pos]['T'] += [tmp_t]

        out    = []
        OUTROW = "{exon},{codon},{base},{pos1},{pos2},{pos3}\n"

        for i in range(0,int(seqlength/3)):
            for base in ['A', 'C', 'G', 'T']:
                out.append(
                    OUTROW.format(
                        exon  = basefile,
                        codon = i,
                        base  = base,
                        pos1  = counter['pos1'][base][i],
                        pos2  = counter['pos2'][base][i],
                        pos3  = counter['pos3'][base][i]
                    )
                )

        return out 
    
    def run(self):

        with Pool(processes= self.threads) as p:

            pre_h = [*p.map(self.horizontal_composition, self.alns )]
            pre_v = [*p.map(self.vertical_composition  , self.alns )]

            # write them down
            with open( self.horcompname, 'w' ) as f:

                f.write("exon,sequence,base,pos1,pos2,pos3\n")
                for ph in pre_h:
                    f.writelines(ph)

            with open( self.vercompname, 'w' ) as f:
                f.write("exon,codon,base,pos1,pos2,pos3\n")
                for pv in pre_v:
                    f.writelines(pv)

# Codonm(
#     threads= 1,
#     alns= ["data/E1718.NT_aligned.fasta_trimmed"]
# ).run()