#!/usr/bin/env python3

import os
import sys
import dendropy
from multiprocessing import Pool
from fishlifeqc.utils import fas_to_dic


class BLCorrelations:
    def __init__(self,
                species_tree_file = None,
                sequences = None,
                suffix = '_constr.tree',
                threads = 1
                ):
        
        self.species_tree_file = species_tree_file
        self.sequences = sequences
        self.suffix = suffix
        
        self.threads = threads

    @property
    def _spps_tree(self):
        """
        Read species tree file 
        from the initial variable
        `self.species_tree_file`

        Returns
        -------
        Tree : str
            Tree str
        """
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

        return tree.as_string(schema = 'newick')

    def _get_taxa(self, sequence):
        sym_c = lambda mystr: mystr.replace(">", "")
        return list(map(sym_c,list(fas_to_dic(sequence))) )

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
        bmfile = os.path.basename( sequence )

        sys.stdout.write("Processing: %s\n" % bmfile)
        sys.stdout.flush()

        tree = dendropy.Tree.get_from_string(
                self._spps_tree, 
                schema = 'newick'
            )

        out_name = sequence + self.suffix
        taxa     = self._get_taxa(sequence)

        tree.retain_taxa_with_labels( taxa )

        tree.write(
            path = out_name,
            schema = 'newick',
            suppress_edge_lengths = True,
            suppress_internal_node_labels = True,
            unquoted_underscores = True
        )

    def get_constraints(self):
        """
        Run in parallel `self._prunning()`
        """
        self._spps_tree
        with Pool(processes = self.threads) as p:
            for seq_file in self.sequences:
                p.apply_async(self._prunning, (seq_file,)).get()

# sequence = ['/Users/ulises/Desktop/GOL/software/fishlifeqc/fishlifeqc/test_tree/E0699.listd_allsets.NT_aligned.fasta_trimmed_renamed']
# file_tree = "/Users/ulises/Desktop/GOL/software/fishlifeqc/prota_all_trimm_noT.ML_spp.tree"
# self = BLCorrelations(species_tree_file=file_tree, sequences=sequence)

