#!/usr/bin/env python3

import os
import sys
import csv
import argparse
import dendropy
from multiprocessing import Pool


class TreeExplore:
    """
    tree manipulator
    """
    def __init__(self, 
                treefiles = None, 
                schema = 'newick',
                collpasebylen = True,
                minlen = 0,
                outfilename = 't_like.csv',
                threads = 1,
                ):
        """
        initial params
        """
        self.treefiles = treefiles
        self.schema = schema
        self.collapsebylen = collpasebylen
        self.minlen = minlen

        self.outfilename = outfilename
        self.threads = threads

    @property
    def _trees(self):
        print('reading trees')
        mytrees = []

        for i in self.treefiles:
            mytrees.append(
                (i, dendropy.Tree.get(path = i, schema = self.schema))
            )

        return mytrees

    def collapse_byLength(self, tree, minlen):
        """
        collapse internal branch length
        by using a minimal length value
        """

        for nd in tree.postorder_edge_iter():
            if nd.length is None:
                continue
            if nd.is_internal() and nd.length <= minlen:
                    nd.collapse()
    
    def collapse_bySupport(self, tree, minsupp):
        """
        collapse internal branch length
        by using a minimal support value
        """

        pass

    def terminal_clades(self, tree):
        """
        gets small values 
        """

        clades = []
        for i in tree:
            cn = i._child_nodes
            if cn:
                t_clade = {}
                for ii in cn:
                    if ii.taxon:
                        t_clade[ii.taxon.label] = round(ii.edge_length, 6)
                if t_clade:
                    clades.append(t_clade)
        return clades

    def _find_Tlikes(self, file_tree):

        bmfile = os.path.basename( file_tree )

        sys.stderr.write("Processing: %s\n" % bmfile)
        sys.stderr.flush()

        tree   = dendropy.Tree.get( path = file_tree, schema = self.schema )

        if self.collapsebylen:
            self.collapse_byLength(tree, self.minlen)

        # print(tree.as_ascii_plot())
        newclades = self.terminal_clades(tree)

        clade_name = []
        for nc in newclades:
            """
            two float transactions
            is the maximun that 
            a float can be bundle
            """
            tmp = []            
            for k,v in nc.items():
                # here v is on the 
                # complete form of 
                # the float
                if v <= self.minlen:
                    tmp += [k]

            if len(tmp) > 1:
                clade_name.append([
                    bmfile, ",".join(tmp)
                ])
                
        return clade_name

    def find_Tlikes(self):
        """
        find T-like terminations 
        in a phylogenetic tree
        """
        with Pool(processes = self.threads) as p:

            preout = []
            
            for file_tree in self.treefiles:
                result = p.apply_async(self._find_Tlikes, (file_tree,))
                preout.append(result)

            # print('processing trees')

            out = [ ['file', 'clade'] ]
            for p in preout:
                if p.get():
                    out.extend( p.get() )

        if len(out) > 1:
            with open(self.outfilename, 'w') as f:
                writer = csv.writer(f)
                writer.writerows(out)


# def getOpts():

#     parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
#                                      description="""

#                         t_like.py [tree files]
                        
#                                      """)
#     parser.add_argument('treefiles',
#                         nargs="+",
#                         help='Filenames')
#     parser.add_argument('-s','--schema',
#                         metavar="",
#                         type= str,
#                         default= "newick",
#                         help='[Optional] Tree format [Default: newick]') 
#     parser.add_argument('-u','--no_collapse',
#                         action="store_false",
#                         help='''[Optional] If selected, no collapse internal branches by length''')
#     parser.add_argument('-m', '--min_len',
#                         metavar = "",
#                         type    = float,
#                         default = 0.000001,
#                         help    = '[Optional] minimun length to collapse internal branch [Default: 0.000001]')
#     parser.add_argument('-o','--outfile',
#                         metavar="",
#                         type= str,
#                         default= "t_like.csv",
#                         help='[Optional] Out filename [Default: t_like.csv]')
#     parser.add_argument('-n', '--threads',
#                         metavar = "",
#                         type    = int,
#                         default = 1,
#                         help    = '[Optional] number of cpus [Default: 1]')                        
#     args = parser.parse_args()
#     return args

# def main():
#     args = getOpts()
#     TreeExplore(
#         treefiles=args.treefiles,
#         schema= args.schema,
#         collpasebylen= args.no_collapse, # default: true
#         minlen=args.min_len,
#         threads= args.threads
#     ).find_Tlikes()

# if __name__ == "__main__":
#     main()

