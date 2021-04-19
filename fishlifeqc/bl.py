#!/usr/bin/env python3

import os
import sys
import csv
import glob
import copy
# import runshell
import dendropy
from multiprocessing import Pool
from fishlifeseq import headers
from fishlifeqc.utils import fas_to_dic, remove_files, runshell, codon_partitions

# import inspect
# import pprint
# def black_box(weird_obj):
#     pprint.pprint(
#         inspect.getmembers( weird_obj, lambda a:not(inspect.isroutine(a)) ),
#         indent= 4
#     )


myos = sys.platform

if myos == 'darwin':
    RAXML = 'raxmlHPC-PTHREADS-SSE3_Darwin_64bit'

elif myos == 'linux' or  myos == 'linux2':
    RAXML = 'raxmlHPC-PTHREADS-SSE3_Linux_64bit'

elif myos == 'win32':
    RAXML = 'raxmlHPC-PTHREADS-SSE3.exe'


class BLCorrelations:

    def __init__(self,
                species_tree_file = None,
                sequences = None,
                suffix  = '_constr.tree',
                prefix = "BL_",
                with_bl = False,
                # cons_trees = None,
                ratio_threshold = 5,
                pearson_threshold = 0.5,
                threads = 1,
                raxml_exe = RAXML,
                evomodel = 'GTRGAMMA',
                codon_partition = True,
                iterations = 2):
        
        self.species_tree_file = species_tree_file
        self.sequences = sequences
        self.suffix = suffix
        self.prefix = prefix
        self.with_bl = with_bl
        # self.cons_trees = cons_trees
        self.ratio_threshold = ratio_threshold
        self.pearson_threshold = pearson_threshold
        
        self.threads = threads

        self.raxml_exe = raxml_exe
        self.evomodel = evomodel
        self.codon_partition = codon_partition
        self.iterations = iterations

    @property
    def __spps_tree__(self):

        try:
            tree = dendropy.Tree.get(
                        path   = self.species_tree_file,
                        schema = 'newick',
                        preserve_underscores = True
                    )

        except dendropy.dataio.newickreader.NewickReader.NewickReaderMalformedStatementError:
            sys.stderr.write("Error reading: %s\n" % self.species_tree_file)
            sys.stderr.flush()
            exit()

        return tree

    def _get_taxa(self, sequence):
        sym_c = lambda mystr: mystr.replace(">", "")
        return list(map(sym_c, headers(sequence)))

    def _prunning(self, sequence):
        """
        Make a backbone from species
        with only taxa found in a fasta files

        Parameters
        ----------
        sequence : dict
            Full file path to source of data.

        Returns
        -------
        None : None
            It writes a newick tree file
            taxa specified in a fasta file
        """
        # sequence = sequence[0]
        tree = copy.deepcopy(self.__spps_tree__)

        out_name = sequence + self.suffix
        taxa     = self._get_taxa(sequence)

        tree.retain_taxa_with_labels( taxa )

        tree.write(
            path = out_name,
            schema = 'newick',
            suppress_edge_lengths = not self.with_bl,
            suppress_internal_node_labels = True,
            unquoted_underscores = True
        )

        return (sequence, out_name)

    def get_constraints(self, quiet = False):
        """
        Run in parallel `self._prunning()`
        """
        # self._spps_tree
        out = []
        with Pool(processes = self.threads) as p:
            for seq_file in self.sequences:

                if not quiet:
                    bmfile = os.path.basename( seq_file )
                    sys.stdout.write("Pruning tree: %s\n" % bmfile)
                    sys.stdout.flush()

                out.append(
                    p.apply_async(self._prunning, (seq_file,)).get()
                )
        return out

    def uRF(self, tree1, tree2):
        """
        Unweighted Robison-Foulds distances

        it might be more efficient
        by not converting dendropy classes
        into strings, but it might need to
        share same TaxonNamespace.

        This method let to compare
        any couple of tree with different
        TaxonNamespace
        """
        tns = dendropy.TaxonNamespace()

        a = dendropy.Tree.get_from_string(
                        tree1.as_string(schema = 'newick'), 
                        schema = 'newick',
                        taxon_namespace=tns)

        b = dendropy.Tree.get_from_string(
                        tree2.as_string(schema = 'newick'), 
                        schema = 'newick',
                        taxon_namespace=tns)

        return (dendropy
                    .calculate
                    .treecompare
                    .symmetric_difference(a,b))

    def root_distance(self, tree):
        taxa  = [i.taxon.label for i in tree.leaf_node_iter()] 
        dists = tree.calc_node_root_distances()
        return { a:b for a,b in zip(taxa, dists) }

    def root_distance_ratios(self, dist_ref, dist_cons):
        mykeys = list(dist_ref)
        out = {}
        for mk in mykeys:
            ref_bl = dist_ref[mk]
            cons_bl = dist_cons[mk]
            out[mk] = round(cons_bl/ref_bl, 2)
        return out

    def filter_dict(self, mdict, threshold = 5):
        return {k:v for k,v in mdict.items() if v >= threshold}

    def terminal_bls(self, tree):
        """
        terminal branch lengths
        """
        df = {}
        for nd in tree.postorder_edge_iter():
            if nd.length is None:
                continue
            taxn = nd.head_node.taxon
            if taxn:
                df[taxn.label] = nd.length

        return df

    def root_ratios(self, ref_tree, cons_tree, threshold):

        # ref_tree, cons_tree, threshold = cp_spps_tree, _cons_tree, self.ratio_threshold 
        # dist_ref  = self.root_distance(ref_tree)
        # dist_cons = self.root_distance(cons_tree)
        dist_ref  = self.terminal_bls(ref_tree) # {A:float(), B:float()}
        dist_cons = self.terminal_bls(cons_tree)
        ratios    = self.root_distance_ratios(dist_ref, dist_cons)

        return self.filter_dict(ratios, threshold=threshold)

    def waiting_times(self, tree):
        """
        branch lengths from tip to 
        root nodes for each taxon
        """
        out = {}
        for nd in tree.preorder_node_iter():
            if nd.is_leaf():
                tmp_label = nd.taxon.label
                tmp_bls   = []
                while True:
                    tmp_len = nd.edge.length

                    if tmp_len is None:
                        break

                    tmp_bls.append(tmp_len)
                    nd = nd._parent_node

                    if nd is None:
                        break

                out[tmp_label] = tmp_bls
        return out

    def align_values(self, table1, table2):
        mykeys = list(table1)
        out = []
        # print(type(table1))
        for mk in mykeys:
            tmp_values = zip(table1[mk], table2[mk])
            for tmp_pair in tmp_values:
                if not tmp_pair in out:
                    out += [tmp_pair]
        return out

    def r_squared(self, x_values, y_values, y_mean, slope, coeff):

        n = len(x_values)
        rss = 0
        tss = 0
        for i in range(n):
            y_pred = coeff + (slope * x_values[i])
            rss += ( y_values[i] - y_pred ) ** 2
            tss += ( y_values[i] - y_mean ) ** 2

        return 1 - (rss/tss)

    def coeff_deter(self, values):

        x = [ i[0] for i in values ]
        y = [ i[1] for i in values ]

        x_mean = sum(x)/len(x)
        y_mean = sum(y)/len(y)

        n = len(x)
        num  = 0
        deno = 0
        for i in range(n):
            num += ( x[i] - x_mean ) * (y[i] - y_mean )
            deno += ( x[i] - x_mean ) ** 2
        # slope
        m = num/deno
        # coeff
        b = y_mean - (m*x_mean)
        # r-squared
        r2 = self.r_squared(x, y, y_mean, m, b)

        return r2 if r2 < self.r2_threshold else None

    def pearson_corr(self, values):

        x = [ i[0] for i in values ]
        y = [ i[1] for i in values ]

        x_mean = sum(x)/len(x)
        y_mean = sum(y)/len(y)

        n = len(x)
        num   = 0
        y_var = 0
        x_var = 0
        for i in range(n):
            num   += ( x[i] - x_mean ) * (y[i] - y_mean )
            y_var += ( y[i] - y_mean ) ** 2
            x_var += ( x[i] - x_mean ) ** 2

        pearson_r = round( num/(( y_var * x_var ) ** 0.5), 7)

        return pearson_r if pearson_r < self.pearson_threshold else None

    def _BL(self, seq_cons_tree_f):

        myseq,_cons_tree_f = seq_cons_tree_f

        bmfile = os.path.basename(_cons_tree_f)
        sys.stdout.write("Processing BL: %s\n" % bmfile)
        sys.stdout.flush()

        cp_spps_tree = copy.deepcopy(self.__spps_tree__)

        _cons_tree   = dendropy.Tree.get(
                        path   = _cons_tree_f, 
                        schema = 'newick',
                        preserve_underscores = True)

        c_taxa = [i.taxon.label for i in  _cons_tree.leaf_node_iter()]
        cp_spps_tree.retain_taxa_with_labels( c_taxa )

        # select the first taxa to 
        # root both trees
        _spps_tree_taxa = [i.taxon.label for i in cp_spps_tree.preorder_node_iter() if i.is_leaf()]
        rooter = _spps_tree_taxa[0]

        # ref tree rooted
        root_ref = cp_spps_tree.mrca( taxon_labels = [rooter] ).edge
        cp_spps_tree.reroot_at_edge(root_ref, length1 = 0, length2 = root_ref.length)

        # constrained tree rooted
        root_cons = _cons_tree.mrca( taxon_labels = [rooter] ).edge
        _cons_tree.reroot_at_edge(root_cons, length1 = 0, length2 = root_cons.length )

        robfould = self.uRF(cp_spps_tree, _cons_tree)

        if robfould:
            sys.stderr.write("Topologies do not match: %s\n" % bmfile)
            sys.stderr.flush()
            # return (False, None, bmfile, None, None)
            return (False, myseq, None, None)

        ref_bls  = self.waiting_times(cp_spps_tree)
        cons_bls = self.waiting_times(_cons_tree)
        bl_table = self.align_values(cons_bls, ref_bls)

        pearson = self.pearson_corr( bl_table )
        ratio   = self.root_ratios(cp_spps_tree, _cons_tree, threshold = self.ratio_threshold)

        remove_files([ _cons_tree_f ])
        # return (True, myseq, bmfile, ratio, pearson)
        return (True, myseq, ratio, pearson)

    def write_down_seq(self, myseq, deletion_list):

        out_name = os.path.basename(myseq) + "_" + self.prefix.replace("_", "")
        aln      = fas_to_dic(myseq)

        with open( out_name, 'w') as f:
            for k,v in aln.items():
                if not k.replace(">", "") in deletion_list:
                    f.write("%s\n%s\n" % (k,v))

    def BrL_reduce(self, tables):

        ratio_table     = [["seq_name", "tip", "ratio"]]
        pearson_r_table = [["seq_name", "pearson_r"]]
        diff_topo       = [['seq_name']]

        # diff_topo.append(['dfas'])

        for table in tables:

            match_topo,myseq,ratio,pearson = table
            bmfile = os.path.basename(myseq)

            if not match_topo:
                diff_topo.append([bmfile])
                continue

            if ratio:
                for k,v in ratio.items():
                    ratio_table.append([ bmfile, k, v ])

            if pearson:
                pearson_r_table.append([ bmfile, pearson ])

            if not pearson:
                deletion_list = ratio if ratio else []
                self.write_down_seq(myseq, deletion_list)

        if len(ratio_table) > 1:
            with open(self.prefix + "ratio.csv", 'w') as f:
                writer = csv.writer(f)
                writer.writerows(ratio_table)

        if len(pearson_r_table) > 1:
            with open(self.prefix + "pearson_corr.csv", 'w') as f:
                writer = csv.writer(f)
                writer.writerows(pearson_r_table)

        if len(diff_topo) > 1:
            with open(self.prefix + "fail_topology.csv", 'w') as f:
                writer = csv.writer(f)
                writer.writerows(diff_topo)

    def __remove_raxmlfs__(self, seq, seq_basename):

        remove_files(
            [
                seq + ".reduced", 
                seq + "_constr.tree", 
                seq + ".stdout",
                seq_basename + ".partitions",
                seq_basename + ".partitions.reduced"
             ]  + \
            glob.glob( "RAxML_log."    + seq_basename + "*" ) +\
            glob.glob( "RAxML_result." + seq_basename + "*" ) +\
            glob.glob( "RAxML_info."   + seq_basename + "*" )
        )

    def __iter_raxml__(self, seq_pruned):
        """
        Failed is a list of lists: [[a], [b], [c], ...].
        It makes an easy export at bl.py command
        """

        raxml_failed = []
        cons_trees   = []

        if self.codon_partition:
            running_msg = "Getting constrained RAxML tree with codon partition for:"
        else:
            running_msg = "Getting constrained RAxML tree for:"

        for seq,pruned in seq_pruned:
            seq_basename = os.path.basename(seq)

            sys.stdout.write( "%s %s\n" % (running_msg, seq_basename))
            sys.stdout.flush()

            suffix = seq_basename + "_raxml"
            # pruned = seq + "_constr.tree"

            if self.codon_partition:
                part_out = seq_basename + ".partitions"
                partition_cmd = "-q %s" % part_out
                codon_partitions(file = seq, outname = part_out, nexus = False)

            else:
                partition_cmd = ""

            cmd = """
                {raxml}         \
                    -p 12345    \
                    -g {pruned} \
                    -m {model}  \
                    -s {seq}    \
                    {partitions}\
                    -n {suffix} \
                    -T {threads}\
                    --silent    \
                    -N {runs}""".format(
                        raxml   = self.raxml_exe,
                        pruned  = pruned,
                        model   = self.evomodel,
                        seq     = seq,
                        partitions = partition_cmd,
                        suffix  = suffix,
                        threads = self.threads,
                        runs    = self.iterations).strip()
            # print(cmd)
            runshell( (cmd.split(), seq + ".stdout"), type = "stdout" )

            self.__remove_raxmlfs__(seq, seq_basename)

            is_there_out = os.path.isfile("RAxML_bestTree."+suffix)
            
            if not is_there_out:
                raxml_failed.append([seq])
                continue

            cons_trees.append( (seq, "RAxML_bestTree." + seq_basename + "_raxml") )

        return (raxml_failed, cons_trees)

    def __get_constr_tree__(self):
        seq_pruned   = self.get_constraints()
        
        sys.stdout.write("\n")
        sys.stdout.flush()

        return self.__iter_raxml__(seq_pruned)

    def BrLengths(self):
        raxml_failed, cons_trees = self.__get_constr_tree__()

        sys.stdout.write("\n")
        sys.stdout.flush()

        if raxml_failed:
            with open( self.prefix + "raxml_failures.txt", 'w') as f:
                writer = csv.writer(f)
                writer.writerows([['seq_name']] + raxml_failed)

        if not cons_trees:
            return None

        with Pool(processes = self.threads) as p:

            preout = []
            for file_tree in cons_trees:
                # result = self._BL(file_tree)
                result = p.apply_async(self._BL, (file_tree,))
                preout.append(result)

            out = [i.get() for i in preout]
            self.BrL_reduce(out)

# _spps_tree_f = "/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/test_tree/prota_all_trimm_noT.ML_spp.tree"
# _cons_tree_f = "/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/test_tree/example/ATP6_model_fixed.tree"
# _cons_tree_f2 = "/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/test_tree/example/E0699_model_inferred_all.tree"
# _cons_tree_f2 = "/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/E0699_model_inferred_all.tree"
# _cons_tree_f = "/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/test_tree/example/E0699_model_fixed.tree" 
# sequences = "/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/test_tree/E0699.listd_allsets.NT_aligned.fasta_trimmed_renamed"
# sequences2 = "/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/test_tree/ATP6.listd_allsets.NT_aligned.fasta_trimmed_renamed"
# seq = sequences

# self = BLCorrelations(
#     species_tree_file = _spps_tree_f,
#     sequences         = [sequences,sequences2]
# )
# self.BrLengths()


