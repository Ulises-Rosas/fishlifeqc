#!/usr/bin/env python

import os
import re
import sys
import csv
import collections

import dendropy
from multiprocessing import Pool
from fishlifeqc.missingdata import Missingdata
from fishlifeqc.t_like      import TreeExplore
from fishlifeqc.bl          import BLCorrelations, RAXML
from fishlifeqc.utils       import (fas_to_dic, 
                                    runshell, 
                                    remove_files,
                                    export_fasta)

myos = sys.platform

if not(myos == 'linux' or  myos == 'linux2' or myos == 'darwin'):
    sys.stderr.write( "This command is not available for %s yet." % myos )
    sys.stderr.flush()
    exit()

SEQMT   = "seqmt"
MAKERMT = "makermt"
CONSEL  = "consel"
CATPV   = "catpv"

class Consel:
    def __init__(self, 
                raxml_exe   = RAXML,
                evomodel    = 'GTRGAMMA',
                seqmt_exe   = SEQMT,
                makermt_exe = MAKERMT,
                consel_exe  = CONSEL,
                catpv_exe   = CATPV,
                out_report = "au_tests.csv",
                threads     = 1):


        self.out_report  = out_report

        self.raxml_exe   = raxml_exe
        self.seqmt_exe   = seqmt_exe
        self.makermt_exe = makermt_exe
        self.consel_exe  = consel_exe
        self.catpv_exe   = catpv_exe

        self.threads  = threads
        self.evomodel = evomodel

    @staticmethod
    def remove_undetermined_chars(seq_file: str, outname: str) -> None:
        """
        RAxML does not allow to have site likehoods for sites 
        without data. Sites completely empty are posible when
        codon-based alignment are performed as only one position
        out of three could be recovered from given sequences.
        """
        in_aln  = fas_to_dic(seq_file)
        out_aln = Missingdata.close_gaps(in_aln, is_codon_aware=False)
        export_fasta(aln = out_aln, outname=outname)

    def _site_likehood(self, seq_tree):
        """
        Site likelihoods are estimated without accounting
        for codon partitions because:
        i)  Empty sites precludes estimations
        ii) The signal of the change when comparing site 
            likelihood using different trees should be 
            recovered because all estimations for the same site
            use the same Q matrix. Then, the signal of change 
            for a site should be present by using either Q matrix 1 
            or Q matrix 2, and so on.
            
        :returns: |str| site_llh_out 
        """

        seq,tree_nHypos = seq_tree
        seq_basename    = os.path.basename(seq)

        # sys.stdout.write("Getting site likelihoods for: %s\n" % seq_basename)
        # sys.stdout.flush()

        site_lnl_out_suffix = tree_nHypos + ".sitelh"

        seq_gap_close  = seq_basename + "_close"
        std_err_holder = tree_nHypos + ".stdout"
        info_carrier   = "RAxML_info." + site_lnl_out_suffix
        site_lnl_out   = "RAxML_perSiteLLs." + site_lnl_out_suffix

        # gaps are closed
        self.remove_undetermined_chars(seq, seq_gap_close)

        cmd = """
            {raxml}\
                -f g\
                -s {seq}  \
                -m {model}\
                -z {constr}\
                -n {suffix}\
                -T {threads}""".format(
                    raxml   = self.raxml_exe,
                    model   = self.evomodel,
                    seq     = seq_gap_close,
                    constr  = tree_nHypos,
                    threads = 1,
                    suffix  = site_lnl_out_suffix).strip()

        runshell( (cmd.split(), std_err_holder), type = "stdout")
        remove_files([ 
            std_err_holder, 
            info_carrier, 
            seq_gap_close,
            seq_gap_close + ".reduced",
            tree_nHypos
        ])

        is_there_out = os.path.isfile(site_lnl_out)

        return None if not is_there_out else site_lnl_out

    def _site_likehood_iter(self, seq_tree_list):
        for seq_tree in seq_tree_list:
            self._site_likehood(seq_tree)

    def _is_incongruent(self, table):
        return False if table['1'] == '1' else True

    def is_incongruent(self, consel_out):

        table = {}
        with open(consel_out, 'r') as f:
            for i in f.readlines():
                line = i.strip()
                if line and not re.findall("(rank|reading)", line):
                    columns = re.split("[ ]+", line)
                    table[columns[1] ] = columns[2]

        return self._is_incongruent(table)
        
    def _consel_pipe(self, seq_tree):
        """
        return failed (list|None) and is_incongruent (bool | None)
        """

        seq,_ = seq_tree
        seq_basename = os.path.basename(seq)

        sys.stdout.write("Running CONSEL for: %s\n" % seq_basename)
        sys.stdout.flush()

        siteout = self._site_likehood(seq_tree)

        if not siteout:
            return ( [seq_basename, "site_likehoods"], None)

        in_noExtension = siteout.replace(".sitelh", "")
        to_parse_table = in_noExtension + ".out"
        msg_holder     = in_noExtension + ".ignore"
        pvalue_file    = in_noExtension + ".pv"

        seqmt_cmd   = "%s --puzzle %s" % (self.seqmt_exe, in_noExtension)
        makermt_cmd = "%s %s" % (self.makermt_exe, in_noExtension)
        consel_cmd  = "%s %s" % (self.consel_exe, in_noExtension)
        catpv_cmd   = "%s %s" % (self.catpv_exe, in_noExtension)

        runshell( (seqmt_cmd.split(),   msg_holder), type = "stdout" )
        runshell( (makermt_cmd.split(), msg_holder), type = "stdout" )
        runshell( (consel_cmd.split(),  msg_holder), type = "stdout" )

        is_there_pv = os.path.isfile(pvalue_file)

        if not is_there_pv:
            return ( [seq_basename, "consel"], None)

        runshell( (catpv_cmd.split(), to_parse_table), type = "stdout" )

        consel_out = ( None, self.is_incongruent(to_parse_table) )

        remove_files([
            siteout,
            in_noExtension + ".mt",
            in_noExtension + ".vt",
            in_noExtension + ".rmt",
            in_noExtension + ".pv",
            in_noExtension + ".ci",
            in_noExtension + ".ignore",
            to_parse_table
        ])

        return consel_out

class Monophyly(TreeExplore, BLCorrelations, Consel):

    def __init__(self,
                 path = '.',
                 fasta_extension = ".fasta",
                 tree_extension = ".tree",
                 taxonomyfile = None,

                 recycle_monos = True, # internal default
                 force_all = False,
                 tgroup = None,
                 raxml_exe = RAXML,
                 evomodel = 'GTRGAMMA',
                 codon_partition = True,
                 iterations = 10,
                 schema = 'newick',
                 collapsebylen = False, # for collapse
                 minlen = 0.000001,     # for collapse
                 collpasebysupp = True, # for collapse
                 minsupp = 0,           # for collapse
                 out_report  = "para_au_tests.csv",
                 threads = 1
                 ):

        self.path = path
        self.fast_ext = fasta_extension
        self.tree_ext = tree_extension

        self.schema          = schema
        self.collapsebylen   = collapsebylen 
        self.minlen          = minlen        
        self.collpasebysupp  = collpasebysupp
        self.minsupp         = minsupp       
        self.taxonomyfile    = taxonomyfile
        self.threads         = threads
        self.out_report      = out_report

        self.outgroup = None

        self.seqmt_exe   = SEQMT  
        self.makermt_exe = MAKERMT
        self.consel_exe  = CONSEL 
        self.catpv_exe   = CATPV  


        # user selected group 
        # for testing paraphyly
        self.tgroup = tgroup

        # selected groups
        self.s_groups = []

        # use already proposed monophyletic groups
        # found at input trees
        self.recycle_monos = recycle_monos

        # force monophyly on all
        # paraphyletic group regardless
        # of self.tgroup
        self.full_force= force_all

        # raxml vars
        self.raxml_exe  = raxml_exe
        self.evomodel   = evomodel
        self.codon_partition = codon_partition
        self.iterations = iterations

        # internal variables
        self.seq_tree = []
        

    def _get_taxa(self, obj, is_nd = True):
        tmp_iter = obj.leaf_iter() if is_nd else obj.leaf_node_iter()
        return [i.taxon.label for i in tmp_iter]

    def _iter_taxa_file(self, mytaxa):
        #TODO: change index to zero after development
        index = 0
        mygroups = {}

        for mt in mytaxa:

            if not self._t_like_taxa.__contains__(mt):
                continue

            tmp_group = self._t_like_taxa[mt][index]
            if not mygroups.__contains__(tmp_group):
                mygroups[tmp_group] = [mt]
            else:
                mygroups[tmp_group] += [mt]    
        
        return mygroups
        
    def _get_groups(self, c_tree):
        mytaxa = self._get_taxa(c_tree, is_nd=False)
        return self._iter_taxa_file(mytaxa)

    def _get_status_groups(self, c_tree, mygroups):
        """
        True if monophyletic
        """
        index = 0

        group_status = []

        for group,taxa in mygroups.items():

            nd      = c_tree.mrca(taxon_labels = taxa)
            nd_taxa = self._get_taxa(nd, is_nd=True)
            is_same = self._is_same_group(nd_taxa, group_indx = index)

            group_status.append( (group, is_same) )

        return group_status

    def deeper_relationships(self, tree, mono_taxa):
        """
        **Experimental**

        get deeper relationship
        between monophyletic groups
        """
        mono_nodes = {}
        for i in list(mono_taxa):
            tmp_node = tree.mrca( taxon_labels =  mono_taxa[i] )
            mono_nodes[i] = tmp_node

        taken = []
        while True:
            avai = set(list(mono_nodes)) - set(taken)

            if not avai:
                break

            a = list(avai)[0]
            g_nd =  mono_nodes[a]
            g_par = g_nd._parent_node

            if g_par is None:
                taken += [a]
                continue

            g_sister = set(g_par._child_nodes)-set([g_nd])

            match_group = []
            for k,v in mono_nodes.items():
                if k == a:
                    continue

                if not g_sister:
                    break

                if v in list(g_sister):
                    match_group.append(k)
                    g_sister -= set([v])

            if not g_sister:
                new_g = match_group + [a]
                for td in new_g:
                    del mono_nodes[td]

                f_m_g = "(%s)" % ",".join(new_g)
                mono_nodes[f_m_g] = g_par
                taken += new_g
            else:
                taken += [a]

        return mono_nodes

    def _force_mono(self, para_groups, target_group, mygroups, full = False):
        targettaxa = []
        monotaxa   = []
        no_parenth = []

        for k,v in mygroups.items():

            under_quote = ["'%s'" % i for i in v]

            if len(under_quote) <= 1:
                no_parenth += under_quote
                continue

            tmp_parenth = ["(%s)" % ",".join(under_quote)]

            if k in para_groups:
                # if full, `target_group`
                # is actually empty.
                # This takes all 
                # paraphyletic groups
                if full:
                    targettaxa += tmp_parenth
                    continue

                if k == target_group:
                    targettaxa +=  tmp_parenth
                    continue
            else:
                if self.recycle_monos:
                    monotaxa += tmp_parenth
                    continue

            no_parenth += under_quote
                
        all_sets = targettaxa + monotaxa + no_parenth

        return all_sets

    def _just_taxa(self, file_tree):

        bmfile = os.path.basename( file_tree )
        # sys.stdout.write("Processing: %s\n" % bmfile)
        # sys.stdout.flush()

        try:
            c_tree = dendropy.Tree.get(
                        path   = file_tree, 
                        schema = self.schema,
                        preserve_underscores = True
                    )
        except dendropy.dataio.newickreader.NewickReader.NewickReaderMalformedStatementError:
            sys.stderr.write("Error reading: %s\n" % bmfile)
            sys.stderr.flush()

        return self._get_taxa(c_tree, is_nd = False)
    
    def _check(self, seq_tree):
        """
        :returns: |tuple| 
        (seq_file1, constr_str1, original_tree_file1) 
        """
        # seq_tree  = self.seq_tree[0]
        seq,file_tree = seq_tree

        bmfile = os.path.basename( file_tree )
        cons_tree_out = bmfile + "_para_forcedCons.tree"

        sys.stdout.write("Finding paraphyletic groups at: %s\n" % bmfile)
        sys.stdout.flush()

        try:
            c_tree = dendropy.Tree.get(
                        path   = file_tree, 
                        schema = self.schema,
                        preserve_underscores = True,
                        rooting='default-rooted'
                    )
        except dendropy.dataio.newickreader.NewickReader.NewickReaderMalformedStatementError:
            sys.stderr.write("Error reading: %s\n" % bmfile)
            sys.stderr.flush()
            return None
                    
        # c_tree = copy.deepcopy(tree)
        if self.collapsebylen or self.collpasebysupp:
            self.collapse(c_tree)

        self._rooting(c_tree)

        mygroups     = self._get_groups(c_tree)
        group_status = self._get_status_groups(c_tree, mygroups) 
        # True if monophyletic
        # [ (spp1 , True ),
        #   (spp2 , False),
        # ]

        para_groups      = [i[0] for i in group_status  if i[1] is False]
        para_groups_taxa = { p:mygroups[p] for p in para_groups  }
        

        if not para_groups:
            return None
        # note: you can have only one 
        # paraphyletic group if this group
        # is being paraphyletic around
        # monophyletic groups
        target_group = ""
        if not self.full_force:
            # if not fully forced,
            # `target_group` must
            # contain a value.
            # Otherwise, `other_group`
            # are all paraphyletic groups
            for g in self.s_groups:
                if g in para_groups:
                    target_group += g
                    break

            if not target_group:
                sys.stderr.write("No groups to test at: %s\n" % bmfile)
                sys.stderr.flush()
                return None

        all_sets = self._force_mono(
            para_groups  = para_groups,
            target_group = target_group,
            mygroups     = mygroups,
            full         = self.full_force
        )

        c_newick = "(%s);" % ",".join(all_sets)
        with open( cons_tree_out, 'w' ) as f:
            f.write(  c_newick + "\n")

        # return (seq, c_newick, bmfile , # + 2)
        return (seq, cons_tree_out, file_tree, target_group, para_groups_taxa)

    def _taxa2groups(self, taxa):

        mytaxa =  set(taxa)
        c_taxa = [(k,len(v)) for k,v in self._iter_taxa_file(mytaxa).items()]
        sorted_taxa =  sorted(c_taxa, key = lambda kv: kv[1], reverse = True)
        return [g for g,_ in sorted_taxa]

    def _readmetadata(self):

        # self.path = "/Users/ulises/Desktop/GOL/software/fishlifeqc/test/para"
        # myrows = []
        # with open(self.metadata, 'r') as f:
        #     reader = csv.reader(f)
        #     for row in reader:
        #         if not row[0].startswith("#"):
        #             myrows.append(row)
        # return myrows
        files = []
        for i in os.listdir(self.path):
            if os.path.isfile( os.path.join(self.path, i ) ):
                files.append(i)

        trees = [i.replace(self.tree_ext , "") for i in files if re.findall(self.tree_ext, i)]
        alns  = [i.replace(self.fast_ext, "") for i in files if re.findall(self.fast_ext, i)]

        base_count = collections.Counter(  trees + alns )
        pairs      = [k for k,v in base_count.items() if v == 2 ]

        if not pairs:
            sys.stderr.write("\nNo alignment-tree name coupling under given file extensions\n\n")
            sys.stderr.write("  alignment extension : %s\n" % self.fast_ext)
            sys.stderr.write("  tree extension      : %s\n" % self.tree_ext)
            sys.stderr.write("\nat '%s' directory\n\n" % self.path)
            sys.stderr.flush()
            exit()

        myrows = []
        for p in pairs:
            aln  = os.path.join(self.path, p + self.fast_ext)
            tree = os.path.join(self.path, p + self.tree_ext)
            myrows.append( [aln, tree] )

        return myrows

    def _raxml_metadata(self, _check_out: list) -> dict:
        """
        Structure:

        {
            'seq_base': {
                'seq_complete' : value 0, 
                'constr_forced': value 1,
                'original_tree': value 2,
                'testing_para' : value 3,
                'all_para'     : value 4
            },
            ...
        }
        """

        metadata = {}
        for line in _check_out:
            my_seq = os.path.basename(line[0])
            metadata[my_seq] = {
                'seq_complete' : line[0], 
                'constr_forced': line[1],
                'original_tree': line[2],
                'testing_para' : line[3],
                'all_para'     : line[4]
            }

        return metadata

    def _two_hypothesis_file(self, cons_trees, raxml_metadata):
        """
        create two hypothesis file

        Lines contain:
        1. Original tree
        2. Full constraint
        """

        for seq,full_constr in cons_trees:

            seqbmfile  = os.path.basename(seq)

            partial_constr = raxml_metadata[seqbmfile]['constr_forced'] # to remove
            original_tree  = raxml_metadata[seqbmfile]['original_tree']

            out_name_two_hypo = seqbmfile + "_Two_Hypothesis"

            raxml_metadata[seqbmfile]['two_hypo_file'] = out_name_two_hypo

            with open(out_name_two_hypo, 'w') as f:
                
                with open(original_tree, 'r') as H1:
                    f.write( H1.readline() )
                
                with open(full_constr, 'r') as H2:
                    f.write( H2.readline() )

            os.remove( full_constr )
            os.remove( partial_constr )

        return raxml_metadata

    def _run_raxml(self, _check_out):
        """
        :returns: |lists| failed, cons_trees

        cons_tree is a list with tuples containing: 
        - Original sequence
        - Two_Hypothesis file (H1 -> original tree,
                               H2 -> constrained)
        ""
        """
        # _check_out = out
        metadata = self._raxml_metadata(_check_out)
        raxml_in = [ ( v['seq_complete'], v['constr_forced'] ) for v in metadata.values() ]

        failed, cons_trees = self.__iter_raxml__(seq_pruned=raxml_in)

        if not cons_trees:
            return (failed, None)
        
        self._two_hypothesis_file(cons_trees, metadata)

        return (failed, metadata)

    def _make_row(self, exon, para_seqs, para, tested):
        
        len_para_seqs = len(para_seqs)

        return list(zip( [exon]   * len_para_seqs  ,
                         para_seqs                 ,
                         [para]   * len_para_seqs  ,
                         [tested] * len_para_seqs  ))

    def _au_test_row(self, metadata, exon):
        # exon = 'exon1'
        all_para      = metadata['all_para']
        all_para_keys = list(all_para)
        tar_para      = metadata['testing_para']

        # exon, sequence, group, tested
        rows = []
        for ap in all_para_keys:

            para_seqs = all_para[ap]
            is_tested = "yes"

            if not self.full_force:
                is_tested = 'yes' if ap == tar_para else 'no'

            rows.extend( self._make_row(exon, para_seqs, ap, is_tested) )

        return rows
                
    def _monophyly_au_test(self, seq_metadata):

        seq,metadata = seq_metadata

        failed,is_incongruent = self._consel_pipe((
                                            metadata['seq_complete' ], 
                                            metadata['two_hypo_file']
                                        ))
        if is_incongruent:
            return ( None, self._au_test_row(metadata, seq) )

        else:
            return ( failed, None )

    def _write_au_table(self, au_table):
        out_file = "au_table.csv" 

        out = [["exon", "sequence", "paraphyletic_group", "au_tested"]]
        out += au_table

        if len(out) > 1:
            with open( out_file , 'w' ) as f:
                writer = csv.writer(f)
                writer.writerows(out)

            sys.stdout.write("\n\nReport of AU tests written at %s\n" % out_file)
            sys.stdout.flush()

    def _write_failures(self, raxml_failures, au_failures):

        out_file = "para_errors.csv" 

        out = [[ "exon", "where" ]]
        if raxml_failures:
            for i in raxml_failures:
                seq_base = os.path.basename(  i[0]  )
                out.append( [ seq_base, 'RAxML' ]  )

        if au_failures:
            out += au_failures

        if len(out) > 1:
            with open( out_file , 'w' ) as f:
                writer = csv.writer(f)
                writer.writerows(out)

            sys.stdout.write("\n\nReport of errors written at %s\n" % out_file)
            sys.stdout.flush()

    def run(self):

        self.seq_tree = self._readmetadata()
        # self.seq_tree = [i for i in self.seq_tree if re.findall("E0072", i[0])]

        with Pool(processes = self.threads) as p:

            if self.full_force:
                self.s_groups = [] if not self.tgroup else [self.tgroup]

            else:
                sys.stderr.write("\nMapping taxa...\r")
                sys.stderr.flush()

                taxa = []
                for _,file_tree in self.seq_tree:
                    preout = p.map_async(self._just_taxa, (file_tree,)).get()[0]
                    # print(preout)
                    if preout:
                        taxa.extend(preout)

                self.s_groups = self._taxa2groups(taxa)
                sys.stderr.write("Mapping taxa...Ok\n\n")
                sys.stderr.flush()

            out = []
            for seq_tree in self.seq_tree:
                preout = p.map_async(self._check, (seq_tree,)).get()[0]
                if preout:
                    out.append(preout)

            sys.stderr.write("\n")
            sys.stderr.flush()

            # pprint.pprint(out)
            rax_failed, rax_metadata = self._run_raxml(_check_out = out)

            au_test_table = []
            au_failed     = []

            if rax_metadata:

                sys.stderr.write("\n")
                sys.stderr.flush()
                
                for seq_metadata in tuple( rax_metadata.items() ):
                    # seq_metadata
                    failed, rows = p.map_async(self._monophyly_au_test, (seq_metadata,)).get()[0]

                    if failed:
                        au_failed.append(failed)

                    if rows:
                        au_test_table.extend(rows)

            self._write_au_table(au_test_table)
            self._write_failures(rax_failed, au_failed)

# taxonomyfile = "/Users/ulises/Desktop/test/mono_tax_file.csv"
# self = Monophyly(
#         path= "/Users/ulises/Desktop/GOL/software/fishlifeqc/test/para",
#         fasta_extension=".fasta",
#         tree_extension=".tree",
#         taxonomyfile = taxonomyfile,

#         collapsebylen = False, # for collapse
#         minlen = 0.000001,     # for collapse
#         collpasebysupp = True, # for collapse
#         minsupp = 0,           # for collapse

#         recycle_monos = True, # internal default
#         force_all = False,    # if True, it will force mono to ALL para
#         tgroup = None,        # target group
#         schema = 'newick',    # internal default
#         threads = 5,        

#         raxml_exe = RAXML,   # internal default
#         evomodel = 'GTRGAMMA',  
#         codon_partition = True,
#         iterations = 1,
# )

# # print(self._readmetadata())