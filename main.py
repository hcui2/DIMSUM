"""
-----------------------------------------------------------------------

Author: Hongzhu Cui (hcui2@wpi.edu)
The main script to run the network propagation and 
subnetwork extraction steps in DIMSUM
 
-----------------------------------------------------------------------
"""
import copy
import operator
import math

import pandas as pd
import networkx as nx

from network_prop import *
from subnetwork_extraction import *
from utils import *


################################################################
#### Processing the interactome 
################################################################

# read in the gene_id/symbol to uniprot id mapping file. 
gene_id_2_uniprot = dict ([line.strip().split('\t')[::-1] for line in open("./data/" + "uniprot_2_gene_id.txt", 'r')])
gene_symbol_2_uniprot = dict ([line.strip().split('\t')[::-1] for line in open("./data/" + "uniprot_2_gene_symbol.txt", 'r')])

# read in interactome data,
G0 = nx.read_edgelist("./data/" + "interactome_edgelist.txt", delimiter = "\t")


################################################################
#### Processing snp-in tool annotation 
################################################################

# read in the annotation file into a dictionary
# also generate a set of mutations.

# mut_effects is a dictionary, the key mutation, 
# the value is a tuple, which is represented as (protein1, protein2, effect)
mut_effects = {}
with open ("./data/" +  "cad_snpin_annotations.txt") as f:
    for line in f:
        prot1, prot2, mut, effect = line.strip().split('\t')

        if mut in mut_effects:
            mut_effects[mut].append((prot1, prot2, effect))
        else:
            mut_effects[mut] = []
            mut_effects[mut].append((prot1, prot2, effect))


################################################################
#### Processing snp file and gene score file
#### Read in gwas dataset
################################################################

# snp_file_name = 'EUR.ASN.CAD.cad.add.160614.website.txt.gz'
# read in snp p-value file into a dictioanry
snp_file_name = "./data/" + "cad_snps.txt.gz"  # replace this 
snp_df = pd.read_csv(snp_file_name, compression='gzip', sep = '\t', header= None, names = ['snp', 'pvalue'])
snp_dict = dict([(snp,float(pvalue)) for snp, pvalue in zip(snp_df.snp, snp_df.pvalue)])

# read in gene score file into a dictionary
gene_file_name = "./data/" + "cad_genescores.txt.gz" 
gene_df = pd.read_csv(gene_file_name, compression='gzip', sep = '\t', header= 0)

gene_id_dict = dict([(gene_id,float(pvalue)) for gene_id, pvalue in zip(gene_df.gene_id, gene_df.pvalue)])
gene_symbol_dict = dict([(gene_symbol,float(pvalue)) for gene_symbol, pvalue in zip(gene_df.gene_symbol, gene_df.pvalue)])

# read in the seed gene list
# seeds_file_name = snp_file_name.replace(".txt.gz", ".SEEDS.txt")
seeds_file_name = "cad_seeds.txt"
full_seeds_file_name = "./data/"  + seeds_file_name
seed_row_list = [line.strip() for line in open(full_seeds_file_name, 'r')]
seed_list = [ line.split('\t')[0] for line in seed_row_list if line.split('\t')[0] != 'n/a']

# translate the gene_id or gene_symbol to uniprot id
prot_dict = {}
for gene_id in gene_symbol_dict:
    if gene_id in gene_id_2_uniprot:
        prot_dict[gene_id_2_uniprot[gene_id]] = gene_symbol_dict[gene_id]

for gene_symbol in gene_id_dict:
     if gene_symbol in gene_symbol_2_uniprot:
        prot_dict[gene_symbol_2_uniprot[gene_symbol]] = gene_id_dict[gene_symbol]

# make -log transformation of the pvalue
node_pvalue_dict = {k: -1*math.log10(v) for k, v in prot_dict.items()}

# make a copy of the original network.
G = copy.deepcopy(G0)

# if the gene not in the interactome remove from the p-value dictionary
to_remove_list = []
for node in node_pvalue_dict.keys():
    if node not in G.nodes():
        to_remove_list.append(node)
for node in to_remove_list:
    del node_pvalue_dict[node]

# if the node in the interactome not in the p-value dictionary
# add a new key, and set the value as 0
for node in G.nodes():
    if node not in node_pvalue_dict:
        node_pvalue_dict[node] = 0

# normalize the p-value, so the node weight will be within [0,1]
max_p = max(node_pvalue_dict.items(), key=operator.itemgetter(1))[1]
min_p = min(node_pvalue_dict.items(), key=operator.itemgetter(1))[1]

normalized_pvalue_dict = {k: (v - min_p)/(max_p - min_p)  for k, v in node_pvalue_dict.items()}

# take a netowrk object and dictionary (gene, p-value)
# annotate the node
for k, v in node_pvalue_dict.items():
    G.node[k]['pvalue'] = v

for k, v in normalized_pvalue_dict.items():
    G.node[k]['norm_p'] = v

# initialize all the edge weight in Graph as 0
nx.set_edge_attributes(G, 0, 'weight0')

mut_count = 0
for mut in mut_effects:
    if mut not in snp_dict:
        continue

    for interaction_effect_combo in mut_effects[mut]:
        prot1, prot2, effect = interaction_effect_combo

        if effect != "Detrimental":
            continue

        # check whether such an interaction exist in the network.
        # if so, update the edge weight.
        if G.has_edge(prot1, prot2) :
            G[prot1][prot2]['weight0'] +=1
        elif G.has_edge(prot2, prot1):
            G[prot2][prot1]['weight0'] += 1
        else:
            continue

for e in G.edges(data = True):
    e[2]['weight'] = sigmoid(e[2]['weight0'])

################################################################
#### network propagation
################################################################

Wprime = normalized_adj_matrix(G, conserve_heat= False, weighted= True)

Y, Ft1 = network_propagation(G, Wprime, num_its=100)

for index, val in Ft1.items():
    update_p = (max_p - min_p) * val + min_p
    G.node[index]['updated_p'] = update_p

################################################################
#### extract subnetowrk
#### output the added genes and the module to local files.
################################################################

run_dimsum(G, seed_list, 100, snp_file_name)


