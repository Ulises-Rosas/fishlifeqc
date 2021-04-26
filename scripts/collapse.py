#!/usr/bin/env python3

import os
import sys
import argparse
import dendropy
from multiprocessing import Pool
from fishlifeqc.t_like import TreeExplore

class Collapse(TreeExplore):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _export_collapse(self, file_tree):

        bmfile = os.path.basename( file_tree )

        sys.stdout.write("Processing: %s\n" % bmfile)
        sys.stdout.flush()

        try:
            tree = dendropy.Tree.get(
                        path   = file_tree, 
                        schema = self.schema,
                        preserve_underscores = True
                    )
        except dendropy.dataio.newickreader.NewickReader.NewickReaderMalformedStatementError:
            sys.stderr.write("Error reading: %s\n" % bmfile)
            sys.stderr.flush()
            return None

        self.collapse(tree)

        tree_str = tree.as_string(schema = 'newick').replace("'", "")

        with open(file_tree + self.suffix, 'w') as f:
            f.write(tree_str)

    def export_collapse(self):

        with Pool(processes = self.threads) as p:
            preout = []
            for file_tree in self.treefiles:
                # self._export_collapse(file_tree)
                results = p.apply_async(self._export_collapse, (file_tree,))
                preout.append(results)

            for i in preout:
                i.get()

def getOpts():

    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""

            Collapse edges in function of support values. 
            

    * Usage:

        $ collapse.py [tree files]

    * Consider also branch lengths to collapse edges:

        $ collapse.py [tree files] -l
    """)
    parser.add_argument('filenames',
                        nargs="+",
                        help='exons names')
    parser.add_argument('-m', '--min_supp',
                        metavar = "",
                        type    = float,
                        default = 0,
                        help    = '[Optional] minimun support value to collapse internal branch [Default: 0]')
    parser.add_argument('-l','--collapse_bylen',
                        action="store_true",
                        help='''[Optional] If selected, collapse internal branches by length''')
    parser.add_argument('-L', '--min_len',
                        metavar = "",
                        type    = float,
                        default = 0.000001,
                        help    = '[Optional] minimun branch length to collapse internal branch [Default: 0.000001]')
    parser.add_argument('-s','--suffix',
                        metavar="",
                        type= str,
                        default= "_collapsed",
                        help='[Optional] Append a suffix at output files [Default: _collapsed]')         
    parser.add_argument('-n', '--threads',
                        metavar = "",
                        type    = int,
                        default = 1,
                        help    = '[Optional] number of cpus [Default: 1]')   

    args = parser.parse_args()
    return args

def main():
    args = getOpts()
    Collapse(
        treefiles = args.filenames, 
        schema = 'newick',
        collpasebylen = args.collapse_bylen,
        minlen = args.min_len,
        collpasebysupp = True,
        minsupp = args.min_supp,
        suffix = args.suffix,
        threads = args.threads
    ).export_collapse()


if __name__ == "__main__":
    main()