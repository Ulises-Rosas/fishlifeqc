
# import glob
import os
import re
import sys
import subprocess
import xml.etree.ElementTree as ET

import matplotlib.pyplot as plt

def runshell(args):
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


def plotcounts(by, outliers, identity, first = 20):

    df = {}
    if by == "exons":
        outname = "mismatchByExonfiles.png"

        for exon,_,_ in outliers:
            if not df.__contains__(exon):
                df[exon] = 1
            else:
                df[exon] += 1

    elif by == "samples":
        outname = "mismatchBySamplefiles.png"

        for _,spps,_ in outliers:
            if not df.__contains__(spps):
                df[spps] = 1
            else:
                df[spps] += 1
    else:
        sys.stderr.write("'%s' is not a valid category\n" % by)
        sys.stderr.flush()
        exit()


    breaks = sorted(df, key=df.get, reverse=False)[0:first]
    vals   = [ df[i] for i in breaks]


    plt.figure(figsize=(10,5))

    plt.barh(breaks, vals,  align='center')
    plt.xlabel('Number of %s' % by)
    plt.title('Mismatch of group at %s%% identity\n first 20' % identity)

    plt.subplots_adjust(left=0.46, right=0.98, top=0.88, bottom=0.11)

    plt.savefig(outname, dpi=300)
    plt.show(block = False)
    plt.close()
