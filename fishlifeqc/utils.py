
# import glob
import os
import re
import subprocess
import xml.etree.ElementTree as ET

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