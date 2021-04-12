
# import glob
import os
import re
import sys
import subprocess
# import matplotlib.pyplot as plt

def runshell(args, type = ""):
    
    if type == "stdout":
        a, b = args
        with open(b, "w") as f:
            subprocess.call(a, stdout= f)

    else:
        p = subprocess.Popen(args)
        p.communicate()

def isfasta(file):

    header = []
    with open(file, "r") as myfile:
        for i in myfile:
            if re.findall("^>", i):
                header += [i]
                break
    if header:
        return True
    else:
        return False


def fas_to_dic(file):
    
    file_content = open(file, 'r').readlines()
    seqs_list   = []
    
    for i in file_content:
        line = i.strip()
        if line: # just if file starts empty
            seqs_list.append(line) 
    
    keys = [] 
    values = []    
    i = 0
    while(">" in seqs_list[i]):
        keys.append(seqs_list[i])
        i += 1 
        JustOneValue = []

        while((">" in seqs_list[i]) == False):
            JustOneValue.append(seqs_list[i]) 
            i += 1

            if(i == len(seqs_list)):
                i -= 1
                break

        values.append("".join(JustOneValue).upper().replace(" ", ""))
        
    return dict(zip(keys, values))

def codon_partitions(file, outname = None, nexus = True):
    # file = "../E1357.unaligned.fasta_listd_lily.NT_aligned.fasta_trimmed_rblastd"
    aln = fas_to_dic(file)
    lengths = set([len(v) for _,v in aln.items()])

    if lengths.__len__() > 1:
        sys.stderr.write("'%s' file has sequences with different lengths\n" % file)
        sys.stderr.flush()
        return file

    aln_length = lengths.pop()

    if not outname:
        outname  = file + ".nex" 
        
    if nexus:
        with open(outname, 'w') as outf:
            outf.write("#nexus\n")
            outf.write("begin sets;\n")
            outf.write("\tcharset pos_1 = %s: 1-%s\\3;\n" % (file, aln_length))
            outf.write("\tcharset pos_2 = %s: 2-%s\\3;\n" % (file, aln_length))
            outf.write("\tcharset pos_3 = %s: 3-%s\\3;\n" % (file, aln_length))
            outf.write("end;\n")
    else:
        with open(outname, 'w') as outf:
            outf.write("DNA, p1 = 1-%s\\3\n" % aln_length)
            outf.write("DNA, p2 = 2-%s\\3\n" % aln_length)
            outf.write("DNA, p3 = 3-%s\\3\n" % aln_length)

    return None

def codon_partitions_nexus(file):
    """
    special case of `codon_partition`. It
    is left as is due to:
    i)  compatibility with `codon_partitions`
    ii) since it accepts only one arguments, it can
        be parallelized right away
    """
    # file = "../E1357.unaligned.fasta_listd_lily.NT_aligned.fasta_trimmed_rblastd"
    aln = fas_to_dic(file)
    lengths = set([len(v) for _,v in aln.items()])

    if lengths.__len__() > 1:
        sys.stderr.write("'%s' file has sequences with different lengths\n" % file)
        sys.stderr.flush()
        return file

    aln_length = lengths.pop()
    outname    = file + ".nex"

    with open(outname, 'w') as outf:
        outf.write("#nexus\n")
        outf.write("begin sets;\n")
        outf.write("\tcharset pos_1 = %s: 1-%s\\3;\n" % (file, aln_length))
        outf.write("\tcharset pos_2 = %s: 2-%s\\3;\n" % (file, aln_length))
        outf.write("\tcharset pos_3 = %s: 3-%s\\3;\n" % (file, aln_length))
        outf.write("end;\n")

    return None

def remove_files(files):

    for f in files:
        try:
            os.remove(f)
        except FileNotFoundError:
            pass

def export_fasta(aln: dict, outname: str):
    """
    This function assumes sequences names
    are already formated (i.e., names starting with ">")
    """
    with open(outname, 'w') as f:
        for k,v in aln.items():
            f.write( "%s\n%s\n" % (k,v))

def _seq_length(aln:dict)-> int:
    return len(next(iter(aln.values())))

def compare_alns(aln1: dict, aln2: dict) -> list:
    # aln1 = {'a':'---cat-aa', 'b':'---cataaa', 'c':'---cataaa'}
    # aln2 = {'a':'cat-aa', 'b':'cataaa', 'c':'cataaa'}
    aln1_len = _seq_length(aln1)
    aln2_len = _seq_length(aln2)

    aln1_gap =  sum([ v.count("-") for v in aln1.values() ])
    aln2_gap =  sum([ v.count("-") for v in aln2.values() ])

    aln1_nheaders = len(list(aln1))
    aln2_nheaders = len(list(aln2))

    return [ aln1_len, aln1_gap, aln1_nheaders,
             aln2_len, aln2_gap, aln2_nheaders  ]

# def plotcounts(by, outliers, identity, first = 20):

#     df = {}
#     if by == "exons":
#         outname = "mismatchByExonfiles.png"

#         for exon,_,_ in outliers:
#             if not df.__contains__(exon):
#                 df[exon] = 1
#             else:
#                 df[exon] += 1

#     elif by == "samples":
#         outname = "mismatchBySamplefiles.png"

#         for _,spps,_ in outliers:
#             if not df.__contains__(spps):
#                 df[spps] = 1
#             else:
#                 df[spps] += 1
#     else:
#         sys.stderr.write("'%s' is not a valid category\n" % by)
#         sys.stderr.flush()
#         exit()


#     breaks = sorted(df, key=df.get, reverse=False)[0:first]
#     vals   = [ df[i] for i in breaks]


#     plt.figure(figsize=(10,5))

#     plt.barh(breaks, vals,  align='center')
#     plt.xlabel('Number of %s' % by)
#     plt.title('Mismatch of group at %s%% identity\n first 20' % identity)

#     plt.subplots_adjust(left=0.46, right=0.98, top=0.88, bottom=0.11)

#     plt.savefig(outname, dpi=300)
#     plt.show(block = False)
#     plt.close()
