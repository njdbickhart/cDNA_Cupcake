__author__ = 'etseng@pacb.com'

"""
ex:
?CAACTTATTGTCCC 0   1   0   0
?CACCT?ATTGTCCC 0   1   0   0
?CA?CTTATTGTCCC 0   1   0   0
G?ACTGCGGAACTTT 0   1   0   0
GTAAT?CG?A?TTCT 0   1   0   0
GTAATGCG?AACTTT 0   2   0   0

we have:
1. list of haplotypes (ex: "AAT", "A?T", "ATG") and the weight of each haplotype
2. list of SNPs at each position
3. user input of ploidity N

start with assumption of 2 allele:
how do we find the center
"""
import os, sys
import networkx as nx
from collections import Counter
from phasing.io.VariantPhaser import Haplotypes

def make_haplotype_counts_into_graph(hap_obj, isoform_tally):
    """
    :param hap_obj: Haplotype object
    :param isoform_tally: list of (isoform, dict of haplotype count), ex: {'PB.45.1': {0:10, 1:20}}
    """
    G = nx.DiGraph()
    for tally in isoform_tally.itervalues():
        for hap_index, count in tally.iteritems():
            hap_str = hap_obj.haplotypes[hap_index]
            n = len(hap_str)
            for i in xrange(n-1):
                # skip any non-variant-assigned locations
                if hap_str[i]=='?' or hap_str[i+1]=='?': continue
                node1 = (i, hap_str[i])
                node2 = (i+1, hap_str[i+1])
                if node1 in G and node2 in G[node1]:
                    G[node1][node2]['weight'] += count
                else:
                    G.add_edge(node1, node2, weight=count)
    return G

def enumerate_allele_candidates(G, cur_node, cur_str, other_str):
    if G.out_degree(cur_node) == 0: # reached a sink
        print cur_str
    else:
        es = G[cur_node].items()
        es.sort(key=lambda x: x[1]['weight'], reverse=True)
        # ex: [((2, 'A'), {'weight': 23}), ((2, 'C'), {'weight': 1})]
        for (node, _dict) in es:
            enumerate_allele_candidates(G, node, cur_str+node[1])


def make_haplotype_graph(haplotype_strings, err_sub, max_diff_allowed):
    G = nx.Graph()
    n = len(haplotype_strings)
    for i in xrange(n): G.add_node(i)
    hap_str_len = len(haplotype_strings[0])
    max_diff = max(max_diff_allowed, err_sub * hap_str_len)
    for i1 in xrange(n-1):
        for i2 in xrange(i1+1, n):
            s1 = haplotype_strings[i1]
            s2 = haplotype_strings[i2]
            sim = sum((s1[k]==s2[k] or s1[k]=='?' or s2[k]=='?') for k in xrange(hap_str_len))
            if (hap_str_len-sim) < max_diff:
                G.add_edge(i1, i2, weight=sim)
    return G


def make_haplotype_counts(isoform_tally):
    """
    :param hap_obj: Haplotype object
    :param isoform_tally: list of (isoform, dict of haplotype count), ex: {'PB.45.1': {0:10, 1:20}}
    """
    hap_count = Counter() # haplotype index --> total count
    for tally in isoform_tally.itervalues():
        for hap_index, count in tally.iteritems():
            hap_count[hap_index] += count
    return hap_count

def error_correct_haplotypes(G, hap_count, hap_obj, isoform_tally):
    cliques = [comm for comm in nx.k_clique_communities(G, 2)]
    # adding all orphans as a clique by itself
    nodes_left = set(G.nodes())

    new_hap_obj = Haplotypes(hap_obj.hap_var_positions, hap_obj.ref_at_pos, hap_obj.count_of_vars_by_pos)
    # for each clique, pick the one that has the most counts
    old_to_new_map = {}
    for i,members in enumerate(cliques):
        for hap_index in members: nodes_left.remove(hap_index)
        stuff = [(hap_index, hap_count[hap_index]) for hap_index in members]
        stuff.sort(key=lambda x: x[1], reverse=True)
        # choose the most abundant haplotype that does NOT have a '?'
        chosen_hap_index = stuff[0][0] # init to first one
        for hap_index, count_ignore in stuff:
            hap_str = hap_obj.haplotypes[hap_index]
            if all(s!='?' for s in hap_str):
                chosen_hap_index = hap_index
                break
        new_hap_index, msg = new_hap_obj.match_or_add_haplotype(hap_obj.haplotypes[chosen_hap_index])
        for x in members: old_to_new_map[x] = new_hap_index
    # look through all leftover nodes, add them if doesn't contain '?'
    for hap_index in nodes_left:
        hap_str = hap_obj.haplotypes[hap_index]
        if all(s != '?' for s in hap_str):
            new_hap_index, msg = new_hap_obj.match_or_add_haplotype(hap_obj.haplotypes[hap_index])
            old_to_new_map[hap_index] = new_hap_index

    # now create a new isoform_tally
    new_isoform_tally = {}
    for k,v in isoform_tally.iteritems():
        new_isoform_tally[k] = Counter()
        for old_hap_index, count in v.iteritems():
            if old_hap_index not in old_to_new_map:
                print >> sys.stderr, "Discarding: {0}".format(hap_obj.haplotypes[old_hap_index])
                continue
            new_hap_index = old_to_new_map[old_hap_index]
            new_isoform_tally[k][new_hap_index] += count
    return old_to_new_map, new_hap_obj, new_isoform_tally




