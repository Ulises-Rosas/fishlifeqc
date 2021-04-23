#!/usr/bin/env python3

import argparse
from multiprocessing import Pool
from fishlifeqc.utils import codon_partitions

def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

            Create codon partition in Nexus format

	* Example:

	    $ codon_partition.py [exon files]
            
            note: partition file name will have the same
                  as the entered alignment + ".nex" extension
""",
                                     epilog="")
    parser.add_argument('filenames',
                        nargs="+",
                        help='Filenames')
    parser.add_argument('-r','--raxml_format',
                    action="store_false",
                    help='''[Optional] If selected, raxml format is created''')
    # parser.add_argument('-s','--suffix',
    #                     metavar="",
    #                     type= str,
    #                     default= "",
    #                     help='[Optional] Suffix for codom composition files [Default = ""]') 
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default = 1]')                        
    args = parser.parse_args()
    return args


class Codon_par:

    def __init__(self, filenames = None, threads = 1, nexus = False):

        self.filenames = filenames
        self.threads = threads
        self.nexus = nexus

    def _codon_partition(self, seq):
        return codon_partitions(seq, outname=None, nexus=self.nexus)

    def run(self):

        with Pool(processes = self.threads) as p:

            no_matched = []
            preout     = []
            for seq in self.filenames:
                results = p.apply_async(self._codon_partition, (seq,))
                preout.append(results)

            for pr in preout:
                gotit = pr.get()
                if gotit:
                    no_matched.append(gotit)

            no_partitions = list(filter(None, no_matched))

            if no_partitions:    
                with open('no_partitions.txt', 'w') as f:
                    f.writelines(no_partitions)

def main():
    args = getOpts()

    Codon_par(
        filenames = args.filenames, 
        threads   = args.threads,
        nexus     = args.raxml_format
        ).run()

if __name__ == "__main__":
    main()