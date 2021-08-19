
# Summary statistics currently available

|variable       |description                                                |note                                                      |
|:--------------|:----------------------------------------------------------|:---------------------------------------------------------|
|aln_base       |alignment name                                             |                                                          |
|nheaders       |number of sequences at the alignment                       |                                                          |
|pis            |parsimony-informative sites                                |                                                          |
|vars           |variable sites                                             |                                                          |
|seq_len        |alignment length                                           |                                                          |
|seq_len_nogap  |number of sites without gaps                               |                                                          |
|gap_prop       |proportion of gaps characters for the alignment matrix     |Considered gap characters: 'N', '-', '!', '?'             |
|nogap_prop     |proportion of non-gap characters for the alignment matrix  |                                                          |
|gc_mean        |mean of GC-content per sequence                            |                                                          |
|gc_var         |variance of GC-content per sequence                        |                                                          |
|gap_mean       |mean of gap percentage per sequence                        |                                                          |
|gap_var        |variance of gap percentage per sequence                    |                                                          |
|pi_mean        |mean of pairwise identity                                  |                                                          |
|pi_std         |standard deviation of pairwise identity                    |                                                          |
|total_tree_len |sum of all branch lengths in tree                          |                                                          |
|treeness       |sum of internal branch lengths devided by 'total_tree_len' |Phillips & Penny 2003, DOI: 10.1016/s1055-7903(03)00057-5 |
|coeffVar_len   |Coefficient of variation of branch lengths after mid-point rooting| Vankan et al. 2021, DOI: 10.1093/sysbio/syab051   |
|inter_len_mean |mean of internal branch lengths in tree                    |                                                          |
|inter_len_var  |variance of internal branch lengths in tree                |                                                          |
|ter_len_mean   |mean of terminal branch lengths in tree                    |                                                          |
|ter_len_var    |variance of terminal branch lengths in tree                |                                                          |
|supp_mean      |mean of support values at nodes in tree                    |                                                          |
|rcv            |relative composition variability                           |Phillips & Penny 2003, DOI: 10.1016/s1055-7903(03)00057-5 |
|treeness_o_rcv |treeness' devided by 'rcv'                                 |Phillips & Penny 2003, DOI: 10.1016/s1055-7903(03)00057-5 |
|saturation     |slope of sequence pairs substitution-distance correlation  |Philippe et al. 2011, DOI: 10.1371/journal.pbio.1000602   |
|SymPval        |symmetry test p-value of the most divergent sequence pair  |Naser-Khdour et al. 2019, DOI: 10.1093/gbe/evz193         |
|MarPval        |marginal test p-value of the most divergent sequence pair  |Naser-Khdour et al. 2019, DOI: 10.1093/gbe/evz193         |
|IntPval        |internal test p-value of the most divergent sequence pair  |Naser-Khdour et al. 2019, DOI: 10.1093/gbe/evz193         |
|LB_std         |standard deviation of the long branch score per sequence   |Struck et al. 2014, DOI: _10.4137/EBo.s14239              |
|tax_prop       |proportion of taxonomical group presence                   |if a taxonomy file with group sets is given               |
|gc_mean_pos1   |mean of GC-content per sequence at codon position 1        |if codon aware option is selected                         |
|gc_var_pos1    |variance of GC-content per sequence at codon position 1    |if codon aware option is selected                         |
|gap_mean_pos1  |mean of gap percentage per sequence at codon position 1    |if codon aware option is selected                         |
|gap_var_pos1   |variance of GC-content per sequence at codon position 1    |if codon aware option is selected                         |
|gc_mean_pos2   |mean of GC-content per sequence at codon position 2        |if codon aware option is selected                         |
|gc_var_pos2    |variance of GC-content per sequence at codon position 2    |if codon aware option is selected                         |
|gap_mean_pos2  |mean of gap percentage per sequence at codon position 2    |if codon aware option is selected                         |
|gap_var_pos2   |variance of GC-content per sequence at codon position 2    |if codon aware option is selected                         |
|gc_mean_pos3   |mean of GC-content per sequence at codon position 3        |if codon aware option is selected                         |
|gc_var_pos3    |variance of GC-content per sequence at codon position 3    |if codon aware option is selected                         |
|gap_mean_pos3  |mean of gap percentage per sequence at codon position 3    |if codon aware option is selected                         |
|gap_var_pos3   |variance of GC-content per sequence at codon position 3    |if codon aware option is selected                         |
|RF             |robinson-fould distance between tree and a reference tree  |if a reference tree is given                              |
