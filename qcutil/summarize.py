

import os
import csv
import sys
import copy
import dendropy
import fishlifeseq
import collections
from multiprocessing import Pool
from fishlifeqc.utils import fas_to_dic

class Stats:

    def __init__(self,
                fastas = None,
                align_based = False,
                seq_based = False,
                prefix = 'stats',
                threads = 1):

        self.gap_chars = ['N', '-', '!']                
        self.fastas    = fastas
        self.threads   = threads
        self.align_based = align_based
        self.seq_based = seq_based
        self.prefix = prefix

    def _stat_sites_std(self, aln):
        """
        # Stats for alignments

        Returning in order:

            * Phylogenetic informative sites
            * Variable sites
            * Proportion of gaps
            * Number of non-gaps characters
            * Length of the alignment
            * Number of sequences

        """
        nheaders = len(aln)
        seq_len = len(next(iter(aln.values())))
        
        pis_s  = 0
        n_gaps = 0
        var_s  = 0
        for pos in range(seq_len):

            column = []
            for v in aln.values():
                column.append( v[pos])

            data_sum = collections.Counter(column)

            for gc in self.gap_chars:
                if data_sum.__contains__(gc):
                    n_gaps += data_sum[gc]
                    del data_sum[gc]

            uniq_char = len(data_sum)

            if not uniq_char:
                continue

            if uniq_char > 1:
                var_s += 1
                tmp_pi = [i for i in data_sum.values() if i > 1]

                if len(tmp_pi) > 1:
                    pis_s += 1

        AREA = nheaders * seq_len

        return (pis_s, var_s, n_gaps/AREA, AREA - n_gaps, seq_len, nheaders)

    def stat_sites(self, fasta_file):

        aln = fas_to_dic(fasta_file)

        if self.align_based:
            al_name = os.path.basename(fasta_file)

            return [al_name] + list(self._stat_sites_std(aln))
        else:
            return self._stat_sites_std(aln)

    def _write_report(self, out, all_spps, all_seqs):
        """
        # Write report:

        * Parsimony informative sites
        * Variable sites                       
        * Number of sites            
        * Number of unique sequences           
        * Gap occupancy (concatenated)   
        * Nucleotide occupancy (concatenated)
        * Mean locus occupancy                 
        * Sum of gaps per locus                        

        ## Mean of locus occupancy percentage

        \nCalculation logic:
        \nmean = (prop_exon1 + prop_exon2 + ...)/number_exons 
        \nmean = (headers_exon1/total_headers + headers_exon2/total_headers + ... )/number_exons
        \nmean = (headers_exon1 + headers_exon2 + ...)/( number_exons * total_headers )
        \nmean = SUM[headers per exon] / (number_exons * total_headers)

        ## Nucleotide ocuppancy (when concatenated)

        \nCalculation logic:
        \nno = All_available_nucleotides/All_possible_nucleotides
        \nno = ( nuc_exon1 + nuc_exon2 + ... )/( total_headers * len_concatenated )
        \nno = SUM[nuc per exon]/( total_headers * (len_exon1 + len_exon2 + ...) )
        \nno = SUM[nuc per exon]/( total_headers * SUM[len per exon] )
        """
        # total headers
        H = len(all_spps)
        N = len(all_seqs) 
        a = [sum(x) for x in zip(*out)]

        # Mean of locus occupancy percentage
        mlop = a[5] * 100 / ( H * N ) 

        # Mean of locus gap percentage
        mlgp = a[2] * 100 / N


        # Nucleotide ocuppancy (when concatenated)
        no   = a[3] * 100 / ( H * a[4] ) 
        go   = 100 - no 

        table = [
            [ "Parsimony informative sites"          , a[0]          ],
            [ "Variable sites"                       , a[1]          ],
            [ "Number of sites"                      , a[4]          ],
            [ "Gap occupancy (concatenated)"         , "%.2f" % go   ],
            [ "Nucleotide occupancy (concatenated)"  , "%.2f" % no   ],
            [ "Number of unique sequences"           , H             ],
            [ "Mean locus occupancy"                 , "%.2f" % mlop ],
            [ "Mean locus gap percentage"            , "%.2f" % mlgp ]
        ]

        col_1_j = "%-36s"
        col_2_j =  "%{}s".format( max( [ len( "%s" % v ) for _,v in table ] ) )

        str_ft = "{}: {}".format( col_1_j, col_2_j )

        print()
        for c,v in table:            
            print(str_ft % (c,v))
        print()

    def _write_report_aln(self, out, all_spps):

        s_table = sorted( out,
                          key = lambda l: l[6], # number of headers
                          reverse = False )
        H = len(all_spps)

        out_table = []
        for na,pis,vars,gap_p,_,seq_len,nhe in s_table:

            out_table.append([
                na,
                "%.4f" % (nhe/H),
                nhe,
                seq_len,
                "%.2f" % (gap_p * 100),
                pis,
                vars
            ])

        col_table = [["aln",
                     "seqs_per_aln_prop",
                     "seqs_per_aln",
                     "sites",
                     "gap_perc",
                     "pis",
                     "var_sites"
                     ]]
     
        out_colum = col_table + out_table

        outfile = self.prefix + "_perAln.tsv"

        with open(outfile, 'w') as f:
            writer = csv.writer(f, delimiter = "\t")
            writer.writerows(out_colum)

        sys.stdout.write( "\nReport for alignments written at %s\n" % outfile )
        sys.stdout.flush()

    def _stats_headers(self, all_spps):
        
        s_table = sorted(
             collections.Counter(all_spps).items(),
             key = lambda kv: kv[1]
             )

        new_table = [["seq",
                      "alns_per_seq_prop",
                      "alns_per_seq",
                      "cum_sum",
                      "cum_sum_per",
                     ]]
     
        nrow    = len(s_table)
        nfastas = len(self.fastas)

        for n,spps_count in enumerate(s_table):

            spps,count = spps_count
            new_table.append([
                spps.replace(">", ""),
                "%.4f" % (count/nfastas),
                count,
                n + 1,
                "%.2f" % ((n + 1)*100/nrow)
            ])

        outfile = self.prefix + "_perSeq.tsv"

        with open(outfile, 'w') as f:
            writer = csv.writer(f, delimiter = "\t")
            writer.writerows(new_table)

        sys.stdout.write("\nReport for sequences written at %s\n" % outfile)
        sys.stdout.flush()

    def run(self):

        with Pool(processes=self.threads) as p:

            preout = []
            for seq in self.fastas:
                result = p.apply_async( fishlifeseq.headers, (seq,) )
                preout.append(result)

            all_spps = []
            for pr in preout:
                all_spps.extend( pr.get() )

            if self.seq_based:    
                self._stats_headers(all_spps)

                if not self.align_based:
                    return None

            all_spps_uniq = list(set(all_spps))

            preout = []
            for seq in self.fastas:
                result = p.apply_async( self.stat_sites, (seq,) )
                preout.append(result)

            out = []
            for pr in preout:
                out.append(  pr.get() )

            if self.align_based:
                self._write_report_aln(out, all_spps_uniq)

            else:
                self._write_report(out, all_spps_uniq, self.fastas)


# import glob

# aln_glob = "/Users/ulises/Desktop/GOL/software/fishlifeqc/demo/mdata/*.fasta"
# # aln_glob = "/Users/ulises/Desktop/GOL/data/alldatasets/nt_aln/internally_trimmed/*_trimmed"
# alns = glob.glob(aln_glob)
# ref_tree = "/Users/ulises/Desktop/GOL/software/fishlifeqc/demo/bl/prota_all_trimm_noT.ML_spp.tree"

# self = Stats(
#     fastas= alns,
#     align_based=True,
#     threads= 3
# )


class Incongruence:

    def __init__(self, 
                gene_trees = None, 
                reference_tree = None, 
                weighted_rf = False,
                quiet = False,
                threads = 1
                ):

        self.gene_trees = gene_trees
        self.reference_tree = reference_tree
        self.weighted_rf = weighted_rf
        self.quiet = quiet
        self.threads = threads

        if self.reference_tree:
            self.__spps_tree__ = dendropy.Tree.get(
                                    path   = self.reference_tree, 
                                    schema = 'newick',
                                    preserve_underscores = True,
                                    rooting='default-rooted'
                                )

    def uRF(self, tree1, tree2):
        """
        # Robison-Foulds distances

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

        if not self.weighted_rf:

            return (dendropy
                        .calculate
                        .treecompare
                        .symmetric_difference(a,b))
        else:
            return (dendropy
                        .calculate
                        .treecompare
                        .weighted_robinson_foulds_distance(a,b))

    def _first_tip(self, tree):
        """
        First tip after arranging tips in preorder
        """
        out = ''
        for i in tree.preorder_node_iter():
            if i.is_leaf():
                out += i.taxon.label
                break

        return out

    def _estimate_rf(self, _gt_tree_f):
        # _cons_tree_f = gene_trees[0]

        bmfile = os.path.basename(_gt_tree_f)

        if not self.quiet:
            sys.stderr.write("Processing: %s\n" % bmfile)
            sys.stderr.flush()


        cp_spps_tree = copy.deepcopy(self.__spps_tree__)

        _gt_tree = dendropy.Tree.get(
                        path   = _gt_tree_f, 
                        schema = 'newick',
                        preserve_underscores = True,
                        rooting='default-rooted')

        gt_taxa = [i.taxon.label for i in  _gt_tree.leaf_node_iter()]
        cp_spps_tree.retain_taxa_with_labels( gt_taxa )

        # select the first taxa to 
        # root both trees
        rooter = self._first_tip(cp_spps_tree)

        # ref tree rooted
        root_ref = cp_spps_tree.mrca( taxon_labels = [rooter] ).edge
        cp_spps_tree.reroot_at_edge(root_ref, length1 = 0, length2 = root_ref.length)

        # gene tree rooted
        root_cons = _gt_tree.mrca( taxon_labels = [rooter] ).edge
        _gt_tree.reroot_at_edge(root_cons, length1 = 0, length2 = root_cons.length )

        # print(cp_spps_tree.as_ascii_plot())
        # print(cp_spps_tree.as_string(schema='newick'))
        return [ bmfile, self.uRF(cp_spps_tree, _gt_tree) ]

    def _rf_report(self, out_table, outfile):

        out_colum = [["gene_tree", "RF_value"]] + out_table

        with open(outfile, 'w') as f:
            writer = csv.writer(f, delimiter = "\t")
            writer.writerows(out_colum)

        if not self.quiet:
            sys.stdout.write( "\nReport written at %s\n" % outfile )
            sys.stdout.flush()

    def get_rf(self, outfile = None):

        out_table =[]
        with Pool(processes=self.threads) as p: 
    
            preout = []
            for gt in self.gene_trees:
                result = p.apply_async(self._estimate_rf, (gt,))
                preout.append(result)

            for pr in preout:
                out_table.append( pr.get() )

        if not outfile:
            return out_table

        else:
            self._rf_report(out_table, outfile)

    def _time_support(self, ref_tree, gene_tree):
        """
        >>> [(time1, support1),(time2, support2)]
        where support is either 1 or -1 when
        node taxa that time is monophyletic
        in the gene tree (support) or not, respectively
        """

        out = []
        for nd in ref_tree.preorder_node_iter():

            if not nd._parent_node:
                continue

            node_taxa = [ i.taxon.label for i in nd.leaf_iter() ]

            gt_nd      = gene_tree.mrca(taxon_labels = node_taxa)
            gt_nd_taxa = [ i.taxon.label for i in gt_nd.leaf_iter() ]

            out.append(( 
                "%.6f" % nd.root_distance, 
                set(node_taxa) == set(gt_nd_taxa),
                node_taxa
            ))

        return out

    def _support_tt(self, _gt_tree_f):

        if not self.quiet:
            bmfile = os.path.basename(_gt_tree_f)
            sys.stderr.write("Processing: %s\n" % bmfile)
            sys.stderr.flush()

        cp_spps_tree = copy.deepcopy(self.__spps_tree__)
        gene_tree = dendropy.Tree.get(
                        path   = _gt_tree_f, 
                        schema = 'newick',
                        preserve_underscores = True
                    )

        gt_taxa = [i.taxon.label for i in  gene_tree.leaf_node_iter()]
        cp_spps_tree.retain_taxa_with_labels( gt_taxa )

        return self._time_support(cp_spps_tree, gene_tree)

    def _format_table(self, table):
        # table = new_table
        s_table = sorted( table.items(),
                          key = lambda kv: float(kv[0]), 
                          reverse =True )

        out = []
        # cum_sum = 0
        for k,v in s_table:
            # change = v['support'] - v['conflict']
            # cum_sum += change
            out.extend(
                [ [ float(k), "supporting_loci" , v['support'],  ",".join(v['support_taxa'] ) ],
                  [ float(k), "conflicting_loci", v['conflict'], ",".join(v['conflict_taxa']) ] ]
                )

            # out.append([
            #     float(k),
            #     v['support'],
            #     v['conflict'],
            #     cum_sum,
            #     ",".join(v['support_taxa']),
            #      ",".join(v['conflict_taxa'])
            # ])

        return out
        
    def reduce_stt_table(self, out):
        # out = out_table

        new_table = {}
        for time,support,taxa in out:

            if not new_table.__contains__(time):

                new_table[time] = {
                    'support'      : 0, 
                    'conflict'     : 0,
                    'support_taxa' : [],
                    'conflict_taxa': []
                }

            if not support:
                new_table[time]['conflict'] += 1
                new_table[time]['conflict_taxa'] += taxa

            else:
                new_table[time]['support'] += 1
                new_table[time]['support_taxa'] += taxa

        # getting unique taxa per time unit 
        # This approach is faster than getting
        # unique taxa in above loop
        for k,v in new_table.items():
            new_table[k]['support_taxa']  = list( set(v['support_taxa'])  )
            new_table[k]['conflict_taxa'] = list( set(v['conflict_taxa']) )

        return self._format_table(new_table)

    def _itt_report(self, out_table, outfile):

        out_table = self.reduce_stt_table(out_table)

        col_table = [["time",
                      "loci_type",
                      "freq",
                      "taxa"]]

        out_colum = col_table + out_table

        with open(outfile, 'w') as f:
            writer = csv.writer(f, delimiter = "\t")
            writer.writerows(out_colum)

        if not self.quiet:
            sys.stdout.write( "\nReport written at %s\n" % outfile )
            sys.stdout.flush()
  
    def support_tt(self, outfile = None):

        self.__spps_tree__.calc_node_root_distances()

        out_table =[]
        with Pool(processes=self.threads) as p: 
    
            preout = []
            for gt in self.gene_trees:
                result = p.apply_async(self._support_tt, (gt,))
                preout.append(result)

            for pr in preout:
                out_table.extend( pr.get() )

        if not outfile:
            return self.reduce_stt_table(out_table)

        else:
            self._itt_report(out_table, outfile)

# import inspect
# import pprint
# def black_box(weird_obj):
#     pprint.pprint(
#         inspect.getmembers( weird_obj, lambda a:not(inspect.isroutine(a)) ),
#         indent= 4
#     )
# import glob
# gene_tree_glob = "/Users/ulises/Desktop/GOL/software/fishlifeqc/demo/para/*.tree"
# gene_trees = glob.glob(gene_tree_glob)
# ref_tree = "/Users/ulises/Desktop/GOL/software/fishlifeqc/demo/bl/prota_all_trimm_noT.ML_spp.tree"

# self = Incongruence(
#     gene_trees = gene_trees,
#     reference_tree = ref_tree,
#     weighted_rf =  True
# )

# table = self.support_tt()