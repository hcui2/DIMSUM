"""
-----------------------------------------------------------------------

Author: Hongzhu Cui(hcui2@wpi.edu)
Implementation of the subnetwork extraction procedure in DIMSUM.  
 
-----------------------------------------------------------------------
"""


def dimsum(G, S, max_iterats =200 ):
    '''
    This core function of the DIMSUM algorithm

    Inputs:
        - G: a network with node weighted by the GWAS and edge weighte by the snp-in tool annotation
        - S: seed gene list. 
        - max_iterats: max number of iterations in the loop, i.e. the number of genes to extract.
    Outputs: 
        - cluster_nodes: all the genes in the disease module
        - added_nodes: the new added genes after the subnetwork extraction procedure. 
    '''

    # Setting up initial set of nodes in cluster

    all_nodes = set(G.nodes())
    cluster_nodes = set(S)
    not_in_cluster = all_nodes - cluster_nodes
    added_nodes = []

    # main loop
    while len(added_nodes) < max_iterats:
        # get cluster neighbors
        neighbor_list = []
        for node in cluster_nodes:
            neighbor_list.extend([neighbor for neighbor in G.neighbors(node)])

        neighbor_set = not_in_cluster & set(neighbor_list)
        if (len(neighbor_set) == 0):
            break

        max_impact = 0
        node_to_add = None

        for candidate in neighbor_set:

            # if candidate not in neighbor_list:
            #     continue;

            # extract updated pvalue
            updated_p = G.node[candidate]['updated_p']

            # calculate the damage
            damage = 1
            for cand_neighor in G.neighbors(candidate):
                if (cand_neighor in cluster_nodes):
                    damage += G[cand_neighor][candidate]['weight0']

            # calculate the impact score.
            iScore = updated_p * damage

            if (iScore > max_impact):
                max_impact = iScore
                node_to_add = candidate

        # remove node from not_cluster
        not_in_cluster.remove(node_to_add)
        # add it to cluster
        cluster_nodes.add(node_to_add)
        # add it to added_nodes
        added_nodes.append(node_to_add)

    return cluster_nodes, added_nodes


def run_dimsum(G, S, iterations, snp_file_name = None):
    '''
    This main function to run the DIMSUM algorithm

    Inputs:
        - G: a network with node weighted by the GWAS and edge weighte by the snp-in tool annotation
        - S: seed gene list. 
        - iterations: max number of iterations in the loop, i.e. the number of genes to extract.
        - snp_file_name: the snp_file_name, used for properly naming the output files. 

    '''
    
    # prepare output files
    outfile1 = "module_gene_list.txt" 
    outfile2 = "added_gene_list.txt"
    if snp_file_name is not None: 
        outfile1 = snp_file_name.replace(".txt.gz", str(iterations)+"_module_gene_list.txt")
        outfile2 = snp_file_name.replace(".txt.gz", str(iterations)+"_added_gene_list.txt")

    # throwing away the seed genes that are not in the network
    all_genes_in_network = set(G.nodes())
    seed_genes = set(S)
    candidate_genes = seed_genes & all_genes_in_network

    if len(candidate_genes) != len(seed_genes):
        print ("ignoring %s of %s seed genes that are not in the network" %(
            len(seed_genes - all_genes_in_network), len(seed_genes)))

    cluster_nodes, added_nodes = dimsum(G, candidate_genes, max_iterats= iterations)

    with open(outfile1, 'w') as f:
        for item in cluster_nodes:
            f.write("%s\n" % item)

    with open(outfile2, 'w') as f:
        for item in added_nodes:
            f.write("%s\n" % item)

    return 



