#!/usr/bin/env python3

import os
import sys
import csv
import collections
# import argparse
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
                 collpasebylen = False,
                 minlen = 0.000001,
                 collpasebysupp = True,
                 minsupp = 0,
                 outfilename = 't_like.csv',
                 taxnomyfile = None,
                 outgroup = None,
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
        self.outgroup = outgroup
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

    @property
    def _t_like_taxa(self):

        if not self.taxonomyfile:
            sys.stderr.write("fishlifeqc t_like: error: the following arguments are required: -t/--taxonomyfile\n")
            sys.stderr.flush()
            exit()

        myrows = []
        with open(self.taxonomyfile, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if not row[0].startswith("#"):
                    myrows.append(row)

        df = {}
        duplications = []
        for row in myrows:
            spps = row[0]

            if df.__contains__(spps):
                duplications.append(spps)

            df[spps] = { n:v  for n,v in enumerate(row[1:])}

        if duplications:
            sys.stderr.write('\nFollowing species names are repeated:\n\n')
            sys.stderr.flush()
            for d in duplications:
                sys.stderr.write(' - %s\n' % d)
                sys.stderr.flush()
            exit()

        return df

    @property
    def _outgroup(self):

        if not self.outgroup:
            return None

        myrows = []
        with open(self.outgroup, 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if not row[0].startswith("#"):
                    if len(row) == 1:
                        myrows.append(["0"] + row)
                    else:
                        myrows.append(row)
        df = {}
        for row in myrows:
            priority = row[0]

            if not df.__contains__(priority):
                df[priority] = [row[1]]
            else:
                if not row[1] in df[priority]:
                    df[priority] += [row[1]]
        return df

    def collapse(self, tree):
        """
        collapse internal branch length
        by using a both branch lengths and
        support values
        """
        # self.collpasebysupp  = True
        # self.collapsebylen = False

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
                # if ed.head_node.label is None:
                #     black_box(ed.head_node)
                if isboth:
                    if ed.length <= minlen:
                        if ed.head_node.label is None:
                            ed.collapse()
                            # continue
                        elif float(ed.head_node.label) <= minsupp:
                            ed.collapse()

                elif islen:
                    if ed.length <= minlen:
                        ed.collapse()

                elif issupp:
                    if ed.head_node.label is None:
                        ed.collapse()
                        # continue
                    elif float(ed.head_node.label) <= minsupp:
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

    def _is_same_group(self,clade_name, group_indx = 0):
        group = []

        for i in clade_name:
            tmp = self._t_like_taxa[i][group_indx]
            if not tmp in group:
                group.append(tmp)

        if len(group) == 1:
            return False if '' in group else True
        
        else:
            return False
        
    def _t_no_samespps(self, clade_names):
        out = []
        # filter same species
        for i in clade_names:
            if not self._is_same_group(i, group_indx = 0):
                out.append(i)

        return out

    def _format_tlike(self, clade_names, bmfile):

        if not clade_names:
            return None

        # clade_names = out
        out = []
        for original,end,p_tax,p_other,reason in clade_names:
            # original,end,reason = clade_names[1]
            p_tax = [] if not p_tax else p_tax
            to_eliminate = set(original) - set(end)
            out.append([
                 bmfile,
                ",".join(original),
                ",".join(p_tax),
                ",".join(to_eliminate),
                ",".join(p_other),
                reason 
            ])

        return out

    def _t_reserve_stem(self, clade, tax_map):

        target_groups = {}
        for t in clade:
            tg = tax_map[t]
            if not target_groups.__contains__(tg):
                target_groups[tg] = [t]
            else:
                target_groups[tg] += [t]

        return target_groups

    def __t_groups__(self, other_taxa, t_clade, tax_map):

        matched_tclade = {}
        matched_oclade = {}

        if not other_taxa:
            return (matched_oclade, matched_tclade)

        in_group = self._t_reserve_stem(t_clade, tax_map)
        out_group = self._t_reserve_stem(other_taxa, tax_map)


        intersection   = set(in_group.keys()) & set(out_group.keys())

        if not intersection:
            return (matched_oclade,matched_tclade)
        
        matched_oclade = {i:out_group[i] for i in intersection}
        matched_tclade = {i:in_group[i] for i in intersection}

        return (matched_oclade, matched_tclade)

    def _t_includes(self, tmp_nd, t_clade_node):
        return [i for i in tmp_nd.postorder_iter() if t_clade_node == i]

    def _t_node_is_mono(self, tmp_nd, tax_map, skip_taxa = None):

        # tmp_nd, tax_map = tmp_node, tax_map
        if tmp_nd.is_leaf():
            tmp_nd = tmp_nd._parent_node

        if not skip_taxa:
            groups = [tax_map[i.taxon.label] for i in tmp_nd.leaf_iter()]
        else:
            groups = []
            for i in tmp_nd.leaf_iter():
                tmp_tip = i.taxon.label
                if not tmp_tip in skip_taxa:
                    groups.append( tax_map[tmp_tip] )

        if len(groups) == 1:
            return False

        return True if len(set(groups)) == 1 else False

    def __recursive_mono__(self, 
                           tree = None,
                           clade = None, 
                           tax_map = None,
                           skip_taxa = None
                           ):

        # t_clade_node = t_clade_node
        # tree         = tree
        # clade        = species
        # tax_map      = tax_map

        are_mono  = []
        to_report = []
        for taxon in clade:
            # taxon
            tmp_node = tree.mrca(taxon_labels=[taxon])
            _is_mono = self._t_node_is_mono(tmp_node, tax_map, skip_taxa)

            if not _is_mono:
                are_mono.append(False)
                to_report.append(taxon)

            else:
                # if not _include_t:
                # return true even if
                # it includes the t clade node
                # if this is considered mono,
                # it will touch t taxa
                # (i.e. P). If the t taxa is 
                # touched, it can be assessed
                # for monophyly
                are_mono.append(True)
            
        if False in are_mono:

            if True in are_mono:

                return (to_report, True)

            else:
                return (to_report, False)
        else:
            return (None, True)

    def _P_tax(self, 
              tree, 
              other_taxa, 
              t_clade, 
              tax_map, 
              skip_taxa = None,
              is_t_mono = None):
        # t_clade_node = t_clade_node
        # tree         = tree
        # other_taxa   = rest_tree
        # t_clade      = t_clade
        # tax_map      = tax_map
        # skip_taxa    = None

        to_report = []
        p_tax     = []

        oclade, iclade = self.__t_groups__(
                                    other_taxa,
                                    t_clade   ,
                                    tax_map   
                                )
        if not oclade:
            return (to_report, p_tax)

        else:
            p_group = []
            for group, species in oclade.items():
                # group, species
                tmp_report, tmp_mono = self.__recursive_mono__(
                                                tree         = tree,
                                                clade        = species, 
                                                tax_map      = tax_map,
                                                skip_taxa    =  skip_taxa
                                            )
                if tmp_mono:
                    # eliminate P taxa
                    # if there is a case
                    # of monophyly for
                    # other taxa
                    p_group.append(group)
                    if tmp_report:
                        to_report.extend(tmp_report)
                else:
                    # if T taxa is not mono
                    # (rare case) delete/report
                    # the P part. 
                    if not is_t_mono:
                        p_group.append(group)
                    else:
                        # if t is mono and not
                        # other taxon (i.e. only tips),
                        # do not eliminate T taxa
                        if tmp_report:
                            to_report.extend(tmp_report)

            if not p_group:
                 return (to_report, p_tax)
            
            for i in p_group:
                p_tax.extend( iclade[i] )

            return (to_report, p_tax)

    def _filter_PD(self, df, t_int_clade, t_clade):
        """
        filter Patristic distances by
        eliminating distances from the same
        t_clade and seperating
        between sister clades
        and internal clades (inside the t_clade)

        Returns
        -------
        /tuple/
        (in_nn, out_nn)
        """
        nearests_in = {}
        nearests_out = {}
        # since keys inside
        # each dictionary are the same,
        # we just take one dictionary key
        # in order to obtain a variable (`tmp_set`)
        # with the precise keys we target on
        if not isinstance( t_clade, list):
            t_clade = list(t_clade)

        tc = t_clade[0] # but is it the closer
                        # to the node?

        tmp_set = set(df[tc]) - set(t_clade)

        if t_int_clade:
            tmp_set -= set(t_int_clade)

            in_tmp = {}
            for tic in t_int_clade:
                in_tmp[tic] = df[tc][tic]

            nearests_in[tc] = in_tmp

        out_tmp = {}
        for toc in tmp_set:
            out_tmp[toc] = df[tc][toc]

        nearests_out[tc] = out_tmp

        out_min = min(nearests_out[tc].values())
        out_nn  = [k for k,v in nearests_out[tc].items() if v == out_min]

        if not nearests_in:
            return (None, out_nn)

        in_min = min(nearests_in[tc].values())
        in_nn  = [k for k,v in nearests_in[tc].items() if v == in_min]

        return (in_nn, out_nn)

    def _NearestNeighborsGroups(self, tree, t_clade, t_clade_node, tax_map):
        '''
        get nearest neighbors groups
        '''
        t_int_clade = []
        # t_clade     = []
        for cn in t_clade_node._child_nodes:
            # a monotipic lineage is a leaf
            if  cn.is_leaf():
                # this leaf is 
                # larger than 
                # self.minlen
                my_leaf = cn.taxon.label
                if not my_leaf in t_clade:
                    t_int_clade.append(my_leaf)
            else:
                for tic in cn.leaf_iter():
                    t_int_clade.append(tic.taxon.label)

        neighbor = [i.taxon for i in t_clade_node.parent_node.leaf_iter()]
        newtree  = tree.extract_tree_with_taxa(taxa = neighbor)
        df = newtree.phylogenetic_distance_matrix().as_data_table()._data

        in_nn, out_nn = self._filter_PD(df, t_int_clade, t_clade)
        out_g = set([tax_map[i] for i in out_nn])

        if not in_nn:
            (None, out_g)

        in_g = set([tax_map[i] for i in in_nn])
        return (in_g, out_g)

    def _t_drop_prop(self, clade, tax_map, prop = 0.5, todrop = None):
        """
        drop taxa and select taxa by the majority of 
        their groups

        returns: an iterable or None
        """
        # clade = t_clade
        # todrop = ['C1']
        # clade = list(tax_map)
        if todrop:
            clade = set(clade) - set(todrop)
            if not clade:
                return clade

        if prop == 1:
            return clade

        check_prop = lambda kv: kv[0] if (kv[1]/len(clade)) >= prop else None

        groups       = collections.Counter([tax_map[i] for i in clade])
        larger_group = list(filter(None, map(check_prop, groups.items())))

        if larger_group:
            return [i for i in clade if tax_map[i] in larger_group]
        else:
            return clade

    def _t_find_stem(self, t_clade):
        # t_clade = [
        #     'Argentiniformes_Argentinidae_Argentina_sp_38',
        #     'Argentiniformes_Argentinidae_Argentina_sp_37'
        # ]
        tax_r  = sorted( list( self._t_like_taxa[t_clade[0]] ), reverse=True)

        stem_r = []
        for tr in tax_r:
            tmp_range = [self._t_like_taxa[i][tr] for i in t_clade]
            uniq_g = set(tmp_range)
            if len(uniq_g) == 1:
                if uniq_g.pop() != "":
                    stem_r.append(tr)
                else:
                    break
            else:
                break

        if not stem_r:
            # use the maximum range
            # possible
            stem_r = [tax_r[0]]

        tax_map = {}
        for k,tax_ranks in self._t_like_taxa.items():
            tmp_range = [tax_ranks[sr] for sr in stem_r]
            tax_map[k] = "-".join(tmp_range)

        return ("-".join(map(str,stem_r)) ,tax_map )

    def _t_reason(self, p_tax, reason):

        if p_tax:
            return reason + ". P_tax: %s" % ','.join(set(p_tax))
        else:
            return reason

    def _t_other_tax(self, node, counter):

        other_taxa = []

        if not node:
            return other_taxa
        # counter = list(filter(None, counter))
        # other_taxa = []

        for l in node.leaf_nodes():
            tmp_label = l.taxon.label
            if not tmp_label in counter:
                other_taxa.append(tmp_label)

        return other_taxa

    def _is_t_monotypic(self, t_clade_node, skip_taxa):
        all_taxons = []

        for spps in t_clade_node.leaf_iter():
            tmp_spps = spps.taxon.label
            if not tmp_spps in skip_taxa:
                all_taxons.append(tmp_spps)

        return True if len(all_taxons) == 1 else False

    def _check_chimps(self, t_clade_node, tax_map, level, P_tax, P_other = None):

        # P_other, P_tax
        passed_nodes = 0
        offset = 0
        # is_mono      = False
        # [i.taxon.label for i in t_clade_node.leaf_iter()]
        is_t_monotypic = self._is_t_monotypic(
                            t_clade_node = t_clade_node,
                            skip_taxa    = P_tax
                        )

        if is_t_monotypic:
            passed_nodes += 1
            offset += 1
            t_clade_node = t_clade_node._parent_node

        source_node = t_clade_node

        while True:
            tmp_mono = self._t_node_is_mono(
                        tmp_nd    = t_clade_node,
                        tax_map   = tax_map, 
                        skip_taxa = P_tax
                    )

            if not tmp_mono:
                break
            
            source_node   = t_clade_node
            t_clade_node  = t_clade_node._parent_node
            passed_nodes += 1

        if passed_nodes == 0 + offset:
            return (False, None, P_other)

        elif passed_nodes == 1 + offset:
            reason = 'T node is monophyletic at %s group' % level

        elif passed_nodes == 2 + offset:
            reason = 'T parent node is monophyletic at %s group' % level

        else:
            reason = 'T parent node and %s node(s) back are monophyletic at %s group' % (passed_nodes - 2, level)
            P_other = set(P_other) - set([i.taxon.label for i in source_node.leaf_iter()])

        if P_tax:
            reason += ' after eliminating P_tax'

        if P_other:
            reason += '. P taxa beyond sister nodes not considered'

        return (True, reason, P_other)        

    def _t_selection(self, tree, t_clade):

        # t_clade = f_clades[2]

        original       = t_clade
        level,tax_map  = self._t_find_stem(t_clade)

        is_t_mono    = len(set([tax_map[i] for i in t_clade])) == 1
        t_clade_node = tree.mrca(taxon_labels=t_clade)

        t_int_tax = self._t_other_tax(t_clade_node, t_clade)
        t_out_tax = self._t_other_tax(t_clade_node._parent_node, t_clade + t_int_tax)
        rest_tree = self._t_other_tax(tree, t_out_tax + t_clade + t_int_tax)

        P_other, P_tax = self._P_tax(tree       = tree, 
                                     other_taxa = rest_tree, 
                                     t_clade    = t_clade, 
                                     tax_map    = tax_map,
                                     skip_taxa  = None,
                                     is_t_mono  = is_t_mono)

        if P_tax:
            t_clade = set(t_clade) - set(P_tax)
            if not t_clade:
                reason = "T clade full of P taxa at %s group" % level
                return (original, t_clade, P_tax, P_other, reason)

        is_mono, reason, P_other = self._check_chimps(t_clade_node = t_clade_node,
                                                      tax_map      = tax_map, 
                                                      level        = level, 
                                                      P_tax        = P_tax, 
                                                      P_other      = P_other)

        if is_mono:
            return (original, t_clade, P_tax, P_other, reason)


        So_other, So_tax  = self._P_tax(tree       = tree, 
                                        other_taxa = t_out_tax, 
                                        t_clade    = t_clade, 
                                        tax_map    = tax_map,
                                        skip_taxa  = P_tax,
                                        is_t_mono  = is_t_mono)

        if So_other:
            P_other += So_other

        Si_other, Si_tax = self._P_tax(tree        = tree, 
                                        other_taxa = t_int_tax, 
                                        t_clade    = t_clade, 
                                        tax_map    = tax_map,
                                        skip_taxa  = P_tax,
                                        is_t_mono  = is_t_mono)

        if Si_other:
            P_other += Si_other

        if not So_tax and not Si_tax:

            reason =  "T clade has no S taxa with any sister node at %s group" % level
            if So_other:
                reason += ". Outside sister P taxa not considered"
            if Si_other:
                reason += ". Internal sister P taxa not considered"

            return (original, t_clade, P_tax, P_other, reason)

        elif not So_tax and Si_tax:

            reason = "T clade has only S taxa with internal sister nodes at %s group" % level

            if So_other:                
                reason += ". Outside sister P taxa not considered"

            tax  = self._t_drop_prop(t_clade, tax_map, prop = 1, todrop=Si_tax)

            # if 0 < len(tax) < len(t_clade):
            #     reason += ". More frequent at shared group level"

            return (original, tax, P_tax, P_other, reason)

        elif So_tax and not Si_tax:

            reason = "T clade has only S taxa with outside sister nodes at %s group" % level

            if Si_other:
                reason += ". Internal sister P taxa not considered"

            tax   = self._t_drop_prop(t_clade, tax_map, prop = 1, todrop=So_tax)

            # if len(So_tax) < len(tax) < len(t_clade):
            #     reason +=". Less frequent at shared group level"

            return (original, tax, P_tax, P_other, reason)

        elif So_tax and Si_tax:
            # neighbor analyses
            in_nn_g, out_nn_g = self._NearestNeighborsGroups(tree, t_clade, t_clade_node, tax_map)
            intersection      = set(in_nn_g) & set(out_nn_g)

            if intersection:
                matched = []
                for tc in t_clade:
                    if tax_map[tc] in intersection:
                        matched.append(tc)

                reason0 = "match with NN groups intersection: '%s'" % ','.join(intersection)
                reason  = reason0 if matched else 'any ' + reason0
                tax     = matched if matched else []

                return (original, tax, P_tax, P_other, reason)

            else:
                o_matched = []
                i_matched = []

                for tc in t_clade:
                    tmp_g = tax_map[tc]

                    if tmp_g in out_nn_g:
                        o_matched.append(tc)

                    if tmp_g in in_nn_g:
                        i_matched.append(tc)

                if o_matched and i_matched:
                    tax    = o_matched + i_matched
                    reason = "match with both NN group"
                    return (original, tax, P_tax, P_other, reason)

                else:
                    tax = self._t_drop_prop(t_clade, tax_map, prop = 1, todrop=So_tax)
                    reason = "shared group (%s) with sister nodes but both were not NN" % level

                    # if len(So_tax) < len(tax) < len(t_clade):
                    #     reason +=". Less frequent at shared group level"

                    return (original, tax, P_tax, P_other, reason)

    def _specific_rooters(self, tree):
        """
        get specific rooters, if available,
        from the gene tree
        """

        all_taxa = [i.taxon.label for i in tree if i.is_leaf()]
        prio_ranges = sorted(map(int, self._outgroup.keys()))

        specific_root = []
        for pr in prio_ranges:
            tmp_range = str(pr)
            tmp_group = self._outgroup[tmp_range]
            tmp_match = [o for o in tmp_group if o in all_taxa]
            if tmp_match:
                specific_root.extend(tmp_match)
                break

        return specific_root

    def _rooting(self, tree):

        if not self._outgroup:
            # print(tree.as_string(schema='newick'))
            # print(tree.as_ascii_plot(plot_metric='length'))
            tree.reroot_at_midpoint(update_bipartitions = True)

        else:
            specific_root = self._specific_rooters(tree)

            if specific_root:
                root_mrca = tree.mrca(taxon_labels=specific_root)

                if len(specific_root) == 1:
                    root_len = root_mrca.edge.length

                    if root_len:
                        tree.reroot_at_edge(root_mrca.edge, 
                                            length1 = 0, 
                                            length2 = root_len, 
                                            update_bipartitions=True, 
                                            suppress_unifurcations= False)

                    else:
                        tree.reroot_at_node(root_mrca,
                                            update_bipartitions=True, 
                                            suppress_unifurcations= False)

                else:
                    tree.reroot_at_node(root_mrca,
                                        update_bipartitions=True, 
                                        suppress_unifurcations= False)

    def _find_Tlikes(self, file_tree):

        # abs_path = "/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/test_tree"
        # tree_path = "E0699.listd_allsets.NT_aligned.fasta_trimmed_renamed.nex.tree"
        # file_tree = os.path.join(abs_path, tree_path)
        # file_tree = '/Users/ulises/Desktop/demo/tlike/E0004.tree'
        # schema = 'newick'
        # file_tree = self.treefiles[0]

        bmfile = os.path.basename( file_tree )

        sys.stdout.write("Processing: %s\n" % bmfile)
        sys.stdout.flush()

        try:
            tree = dendropy.Tree.get(
                        path   = file_tree, 
                        schema = self.schema,
                        preserve_underscores = True,
                        rooting='default-rooted' # makes figtree-like rooting
                    )
        except dendropy.dataio.newickreader.NewickReader.NewickReaderMalformedStatementError:
            sys.stderr.write("Error reading: %s\n" % bmfile)
            sys.stderr.flush()
            return None

        if self.collapsebylen or self.collpasebysupp:
            # collapse affects rooting
            # then, this goes first
            self.collapse(tree)

        self._rooting(tree)
        # print(tree.as_ascii_plot(plot_metric = 'length') )
        # print(tree.as_ascii_plot() )
        # print(tree.as_string(schema = 'newick') )
        
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
                """
                here v is on the complete 
                form of  the float 
                """
                if v <= self.minlen:
                    tmp += [k]

            if len(tmp) > 1:
                clade_name.append(tmp)

        f_clades = self._t_no_samespps(clade_name)

        out = []
        if f_clades:
            for i in f_clades:
                # print(i)
                out.append(self._t_selection(tree, i))

        return self._format_tlike(out, bmfile)

    def find_Tlikes(self):
        """
        find T-like terminations 
        in a phylogenetic tree
        """
        self._t_like_taxa
        self._outgroup
        # print(self._outgroup)

        with Pool(processes = self.threads) as p:

            preout = []
            for file_tree in self.treefiles:
                # result = self._find_Tlikes(file_tree)
                result = p.apply_async(self._find_Tlikes, (file_tree,))
                preout.append(result)

            out = [ ['file','original','P_tax','delete_inside_t', 'delete_outside_t','reason'] ]
            for p in preout:
                gotit = p.get()
                if gotit:
                    out.extend( gotit )

        if len(out) > 1:
            with open(self.outfilename, 'w') as f:
                writer = csv.writer(f)
                writer.writerows(out)

            sys.stdout.write( "\nReport written at %s\r" % self.outfilename )
            sys.stdout.flush()

            self._expand_tlike_table(out)

            sys.stdout.write("\n")
            sys.stdout.flush()

    def _expand_tlike_table(self, tlike_out):

        dirname       = os.path.dirname(self.outfilename)
        basename      = os.path.basename(self.outfilename)
        expanded_name =  os.path.join(dirname, "collapsed_" + basename)

        myrows = [('file', 'sequences', 'delete_type')]

        for row in tlike_out[1::]:

            file      = row[0]
            del_in    = row[3]
            del_out   = row[4]
            del_clade = []

            if del_in:
                del_clade.extend(
                    [ (file, i, 'delete_inside_t') for i in del_in.split(",") ]
                )
            if del_out:
                del_clade.extend(
                    [ (file, i, 'delete_outside_t') for i in del_out.split(",") ]
                )

            if del_clade:
                myrows.extend(del_clade)


        if len(myrows) > 1:
            with open(expanded_name, 'w') as f:
                writer = csv.writer(f)
                writer.writerows(myrows)

            sys.stdout.write( "Report written at %s and %s\n" % (self.outfilename, expanded_name) )
            sys.stdout.flush()
                    
# self = TreeExplore(
#     taxnomyfile= "/Users/ulises/Desktop/GOL/software/fishlifeqc/demo/tlike/taxa_file.csv", 
#     outgroup= None,
#     collpasebysupp=True,
#     treefiles=['/Users/ulises/Desktop/E0713.listd_allsets.NT_aligned.fasta_trimmed.nex.tree']
# )
# self.find_Tlikes()
