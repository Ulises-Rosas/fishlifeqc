#!/usr/bin/env python3

import os
import glob
import shutil
import argparse


def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""
			Pack (or unpack) files into directories

	* Pack files into 40 directories:

	    $ splitexonfiles.py [exon files] -w . -n 40 

	* Unpack directories:

	    $ splitexonfiles.py -u partition* -w .
	 
	    note: The first command packs into directories
		  with the following pattern: `partition[0-9]+`.
		  That's why in this second command we use
		  `partition*` glob pattern for covering all 
                  previously created directories. 
""",
                                     epilog="")
    parser.add_argument('glob',
                        nargs="+",
                        help='Filenames')
    parser.add_argument('-u','--unpack',
                        action="store_true",
                        help='''[Optional] If selected, unpack files from directories into a specific dir''')
    parser.add_argument('-w','--where',
                        metavar="",
                        type= str,
                        default= ".",
                        help='[Optional] Directory where files are moved [Default = "."]')
    parser.add_argument('-n','--partition',
                        metavar="",
                        type= int,
                        default= 3,
                        help='[Optional] Number of partitions [Default = "3"]')
    parser.add_argument('-p','--prefix',
                        metavar="",
                        type= str,
                        default= "partition",
                        help='[Optional] Prefix for each packing directory [Default = "partition"]') 
    args = parser.parse_args()
    return args

def pack(myfiles, npart, where, prefix):
    # npart   = 3
    # myfiles = glob.glob("../data/*")
    # where   = "../data"
    window  = len(myfiles)/npart if len(myfiles) >= npart else 1
    init    = 0

    for i in range(0, npart):
        names = myfiles[round(init): round(init + window)]

        newdir = os.path.join(where, "%s%s" % (prefix, i))
        os.mkdir(newdir)

        for f in names:
            if not os.path.isdir(f):
                shutil.copy(f, newdir)
                os.remove(f)

        init += window

def unpack(mydirs, where):

    for i in mydirs:

        if os.path.isdir(i):
            infiles = glob.glob(os.path.join(i, "*"))

            for f in infiles:
                shutil.copy(f, where)

            shutil.rmtree(i)

def main():
    args = getOpts()

    if args.unpack:
        unpack(mydirs = args.glob,
               where  = args.where)

    else:
        pack(myfiles = args.glob, 
             npart   = args.partition,
             where   = args.where,
             prefix  = args.prefix)
if __name__ == "__main__":
    main()
