#!/usr/bin/env python3

import os
import sys
import csv
import collections
import argparse
import dendropy
from multiprocessing import Pool



# import inspect
# import pprint
# def black_box(weird_obj):
#     pprint.pprint(
#         inspect.getmembers( weird_obj, lambda a:not(inspect.isroutine(a)) ),
#         indent= 4
#     )


class TreeExplore:
    """
    tree manipulator
    """
    def __init__(self, 
                treefiles = None, 
                schema = 'newick',
                collpasebylen = True,
                minlen = 0,
                collpasebysupp = True,
                minsupp = 0,
                outfilename = 't_like.csv',
                taxnomyfile = None,
                suffix = '_fishlife',
                threads = 1):
        """
        initial params
        """
        self.treefiles = treefiles
        self.schema = schema
        self.collapsebylen = collpasebylen
        self.minlen = minlen

        self.collpasebysupp = collpasebysupp
        self.minsupp = minsupp

        self.outfilename = outfilename
        self.threads = threads

        self.taxonomyfile = taxnomyfile
        self.suffix = suffix

    @property
    def _taxa(self):

        if not self.taxonomyfile:
            return None

        myrows = []
        with open(self.taxonomyfile, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                myrows.append(row)

        newnames = [i[1] for i in myrows]

        duplications = []
        for k,v in collections.Counter(newnames).items():
            if v > 1:
                duplications.append(k)

        if duplications:
            sys.stderr.write('\nFollowing new names are repeated:\n\n')
            sys.stderr.flush()
            for d in duplications:
                sys.stderr.write(' - %s\n' % d)
                sys.stderr.flush()
            exit()

        return { on[0]:on[1] for on in myrows  }

    def collapse(self, tree):
        """
        collapse internal branch length
        by using a both branch lengths and
        support values
        """

        # , bylen, minlen, bysupp, minsupp,
        minsupp = self.minsupp
        minlen  = self.minlen

        isboth  = self.collpasebysupp and self.collapsebylen
        islen   = not self.collpasebysupp and self.collapsebylen
        issupp  = self.collpasebysupp and not self.collapsebylen

        for ed in tree.postorder_edge_iter():

            if ed.length is None:
                continue

            if ed.is_internal():

                if isboth:
                    if ed.head_node.label is None:
                        continue

                    if float(ed.head_node.label) <= minsupp and ed.length <= minlen:
                        ed.collapse()

                elif islen:
                    if ed.length <= minlen:
                        ed.collapse()

                elif issupp:
                    if ed.head_node.label is None:
                        continue

                    if float(ed.head_node.label) <= minsupp:
                        ed.collapse()

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

        # file_tree = "../E0012.listd_allsets.NT_aligned.fasta_trimmed.nex.treefile_renamed"
        # schema = 'newick'

        bmfile = os.path.basename( file_tree )

        sys.stderr.write("Processing: %s\n" % bmfile)
        sys.stderr.flush()

        tree   = dendropy.Tree.get( path = file_tree, schema = self.schema )

        if self.collapsebylen or self.collpasebysupp:
            self.collapse(tree)

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
                # print(k,v)
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

            out = [ ['file', 'clade'] ]
            for p in preout:
                if p.get():
                    out.extend( p.get() )

        if len(out) > 1:
            with open(self.outfilename, 'w') as f:
                writer = csv.writer(f)
                writer.writerows(out)

    def _rename_tips(self, file_tree):
        # file_tree = '../E0012.listd_allsets.NT_aligned.fasta_trimmed.nex.treefile'

        bmfile = os.path.basename( file_tree )
        sys.stderr.write("Processing: %s\n" % bmfile)
        sys.stderr.flush()

        tree = dendropy.Tree.get(
                    path   = file_tree,
                    schema = self.schema,
                    preserve_underscores = True)

        tree.as_string(schema = self.schema)

        name_set = self._taxa

        for i in range(len(tree.taxon_namespace)):
            tipname = tree.taxon_namespace[i].label
            if name_set.__contains__(tipname):
                tree.taxon_namespace[i].label = name_set[tipname]

        tree.write(path = file_tree + self.suffix, schema = 'newick')

    def rename_tips(self):

        if not self.taxonomyfile:
            sys.stderr.write("No taxonomy file provided\n")
            sys.stderr.flush()
            return None

        self._taxa

        with Pool(processes = self.threads) as p:
            for file_tree in self.treefiles:
                p.apply_async(self._rename_tips, (file_tree,)).get()

# self = TreeExplore(taxnomyfile= "./../name_map.csv", treefiles=['./../E0012.listd_allsets.NT_aligned.fasta_trimmed.nex.treefile'])
# def getOpts():

#     parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
#                                      description="""

#                         t_like.py [tree files]
                        
#                                      """)
#     parser.add_argument('treefiles',
#                         nargs="+",
#                         help='Filenames')
#     parser.add_argument('-f','--format',
#                         metavar="",
#                         type= str,
#                         default= "newick",
#                         help='[Optional] Tree format [Default: newick]') 
#     parser.add_argument('-l','--collapse_bylen',
#                         action="store_true",
#                         help='''[Optional] If selected, collapse internal branches by length''')
#     parser.add_argument('-L', '--min_len',
#                         metavar = "",
#                         type    = float,
#                         default = 0.000001,
#                         help    = '[Optional] minimun branch length to collapse internal branch [Default: 0.000001]')
#     parser.add_argument('-s','--collapse_bysupp',
#                         action="store_true",
#                         help='''[Optional] If selected, collapse internal branches by support value''')
#     parser.add_argument('-S', '--min_supp',
#                         metavar = "",
#                         type    = float,
#                         default = 0,
#                         help    = '[Optional] minimun support value to collapse internal branch [Default: 0]')
#     parser.add_argument('-o','--outfile',
#                         metavar="",
#                         type= str,
#                         default= "t_like.csv",
#                         help='[Optional] Out filename [Default: t_like.csv]')
#     parser.add_argument('-t','--taxonomyfile',
#                         metavar="",
#                         type= str,
#                         default= None,
#                         help='[Optional] Taxonomy file [Default: None]') 
#     parser.add_argument('-a','--suffix',
#                         metavar="",
#                         type= str,
#                         default= "_fishlife",
#                         help='[Optional] Append a suffix at output files [Default: _fishlife]')         
#     parser.add_argument('-n', '--threads',
#                         metavar = "",
#                         type    = int,
#                         default = 1,
#                         help    = '[Optional] number of cpus [Default: 1]')          
#     args = parser.parse_args()
#     return args

# def main():
#     args = getOpts()
#     # print(args)
#     TreeExplore(
#         treefiles      = args.treefiles,
#         schema         = args.format,
#         collpasebylen  = args.collapse_bylen, # default: false
#         minlen         = args.min_len,
#         collpasebysupp = args.collapse_bysupp, # default: false
#         minsupp        = args.min_supp,
#         taxnomyfile    = args.taxonomyfile,
#         suffix         = args.suffix, 
#         threads        = args.threads
#     ).rename_tips()

# if __name__ == "__main__":
#     main()

