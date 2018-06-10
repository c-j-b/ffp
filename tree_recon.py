import dendropy
import sys
from optparse import OptionParser
import re
import numpy as np


# not used
# from matrix to a dict of dict where key is row col number
def matrix2dict(M, nsameflag=True):
    if type(M) == np.matrixlib.defmatrix.matrix:
        sp = M.shape
        if len(sp) > 2:
            print("Dimension too high!")
            return None
        else:
            a = sp[0]
            b = sp[1]
        if a == 1:
            dic = dict()
            for idx in range(0, b):
                dic[idx] = M[idx]
        elif a > 1:
            dic = dict((x, dict()) for x in range(0, a))
            for idx in range(0, a):
                for idx2 in range(0, b):
                    if nsameflag:
                        if idx1 == idx2:
                            continue
                    dic[idx][idx2] = M[idx, idx2]
    elif isinstance(M, list):
        a = len(M)
        b = len(M[0])
        if a == 1:
            dic = dict()
            for idx in range(0, b):
                dic[idx] = M[idx]
        elif a > 1:
            dic = dict((x, dict()) for x in range(0, a))
            for idx in range(0, a):
                for idx2 in range(0, b):
                    if nsameflag:
                        if idx == idx2:
                            continue
                    dic[idx][idx2] = M[idx][idx2]
    else:
        print("invalid input!")
        return None
    return dic


# from square matrix to a dict of dict where key is in the nl(namelist)
def sqmatrixandnl2dict(M, nl, nsameflag=True):
    if type(M) == np.matrixlib.defmatrix.matrix:
        sp = M.shape
        if len(sp) > 2:
            print("Dimension too high!")
            return None
        else:
            a = sp[0]
        if a == 1:
            dic = {nl[0]: {nl[0]: M[0, 0]}}
        elif a > 1:
            dic = dict((x, dict()) for x in nl)
            for idx1, n1 in enumerate(nl):
                for idx2, n2 in enumerate(nl):
                    if nsameflag:
                        if idx1 == idx2:
                            continue
                    dic[n1][n2] = M[idx1, idx2]
    elif isinstance(M, list):
        a = len(M)
        if a == 1:
            dic = {nl[0]: {nl[0]: M[0]}}
        elif a > 1:
            dic = dict((x, dict()) for x in nl)
            for idx1, n1 in enumerate(nl):
                for idx2, n2 in enumerate(nl):
                    if nsameflag:
                        if idx1 == idx2:
                            continue
                    dic[n1][n2] = M[idx1][idx2]
    else:
        print("invalid input!")
        return None
    return dic


# find the minimum 2 nodes in Q dict to join
def findmin(Q):
    minvalue = 0
    nl = list(Q.keys())
    for i, k1 in enumerate(nl):
        for k2 in nl[i + 1:]:
            if Q[k1][k2] < minvalue:
                m1 = k1
                m2 = k2
                minvalue = Q[k1][k2]
    return m1, m2


# from distance dict to 2 tree dict: treeparentnodes and treebranchlens
def distancedict2treedict_nj(DistanceDict):
    # NeighborJoining
    treeparentnodes = {}
    treebranchlens = {}
    lf_nodes = list(DistanceDict.keys())
    NUM_lf = len(lf_nodes)
    nodeset = set(lf_nodes)
    for inidx in range(0, NUM_lf - 2):
        print(nodeset)
        SUMDM = {}  # sum of distance matrix in each row
        for k in nodeset:
            SUMDM[k] = 0
            for k2 in nodeset - set([k]):
                SUMDM[k] += DistanceDict[k][k2]
        NUM_nodes = len(nodeset)
        Q = dict((x, dict()) for x in nodeset)
        # calculate Q matrix to find 2 nodes to join
        for n1 in nodeset:
            for n2 in nodeset - set([n1]):
                Q[n1][n2] = (NUM_nodes - 2) * DistanceDict[n1][n2] - SUMDM[n1] - SUMDM[n2]
        print(Q)
        m1, m2 = findmin(Q)
        # distance of 2 joining nodes to new parent
        branchlen1 = DistanceDict[m1][m2] / 2.0 + 1.0 / float(2 * (NUM_nodes - 2)) * (SUMDM[m1] - SUMDM[m2])
        branchlen2 = DistanceDict[m1][m2] - branchlen1

        new = 'In_node' + str(inidx)
        treeparentnodes[m1] = new
        treeparentnodes[m2] = new
        treebranchlens[m1] = branchlen1
        treebranchlens[m2] = branchlen2

        nodeset -= set([m1, m2])
        # update distance of new node to other remaining nodes
        DistanceDict[new] = {}
        for k in nodeset:
            newdist = (DistanceDict[m1][k] + DistanceDict[m2][k] - DistanceDict[m1][m2]) / 2.0
            DistanceDict[new][k] = newdist
            DistanceDict[k][new] = newdist
        nodeset = nodeset | set([new])

    # last 2 joining and root
    nodelist = list(nodeset)
    m1, m2 = nodelist[0], nodelist[1]
    treeparentnodes[m1] = 'root'
    treeparentnodes[m2] = 'root'
    treebranchlens[m1] = DistanceDict[m1][m2] / 2.0
    treebranchlens[m2] = DistanceDict[m1][m2] / 2.0  # how to root?
    treeparentnodes['root'] = ''
    treebranchlens['root'] = 0.0

    return treeparentnodes, treebranchlens


# from 2 tree dict: treeparentnodes and treebranchlens to a dendropy.tree object
def TPNandTBL2Tree(treeparentnodes, treebranchlens, lf_nodeslist):
    taxon_namespace = dendropy.TaxonNamespace(lf_nodeslist)
    tree = dendropy.Tree(taxon_namespace=taxon_namespace)
    In_nodes = []
    namenodedict = {'root': tree.seed_node}
    for i in range(0, len(lf_nodeslist) - 2):
        In_nodes.append(dendropy.Node(edge_length=treebranchlens['In_node' + str(i)]))
        namenodedict['In_node' + str(i)] = In_nodes[i]
    lf_nodes = []
    for i, lf in enumerate(lf_nodeslist):
        lf_nodes.append(dendropy.Node(edge_length=treebranchlens[lf]))
        lf_nodes[i].taxon = taxon_namespace.get_taxon(str(lf_nodeslist[i]))
        namenodedict[lf] = lf_nodes[i]
    # all_childnodes=In_nodes+lf_nodes

    activenodes = ['root']
    while len(activenodes) > 0:
        for an in activenodes:
            clist = [c for c, p in treeparentnodes.items() if p == an]
            c1, c2 = clist[0], clist[1]
            namenodedict[an].add_child(namenodedict[c1])
            namenodedict[an].add_child(namenodedict[c2])
            if c1 in treeparentnodes.values():  # is a in_node
                activenodes.append(c1)
            if c2 in treeparentnodes.values():  # is a in_node
                activenodes.append(c2)
            activenodes.remove(an)

    print(tree.as_string("newick"))
    print(tree.as_ascii_plot())
    return tree


# test matrix to recon the tree and reroot tree at edge number rootat and division factor of that edge alpha
def test(rootat, alpha):
    DistanceMatrix = np.matrix([[0, 5, 9, 9, 8],
                                [5, 0, 10, 10, 9],
                                [9, 10, 0, 8, 7],
                                [9, 10, 8, 0, 3],
                                [8, 9, 7, 3, 0]])
    lf_nodeslist = ['A', 'B', 'C', 'D', 'E']  # list(range(0,5))

    DistanceDict = sqmatrixandnl2dict(DistanceMatrix, lf_nodeslist)
    print(DistanceDict)
    assert (DistanceDict != None)
    tpn, tb = distancedict2treedict_nj(DistanceDict)
    print(tpn)
    print(tb)
    tree = TPNandTBL2Tree(tpn, tb, lf_nodeslist)
    eglst = [eg for eg in tree.postorder_edge_iter()]  # can change the iter method to change order of edges
    rootedge = eglst[rootat]
    rel = rootedge.length
    tree.reroot_at_edge(rootedge, length1=rel * alpha, length2=rel * (1 - alpha))
    print(tree.as_string("newick"))
    print(tree.as_ascii_plot())


# not used
def build_in_yield_tree_from_files():
    taxon_namespace = dendropy.TaxonNamespace()
    f1 = open("path/to/trees1.nex", "r")
    f2 = open("path/to/trees2.nex", "r")
    tree_yielder = dendropy.Tree.yield_from_files(
        files=[f1, f2, "path/to/trees3.nex", "path/to/trees4.nex"],
        schema="nexus",
        taxon_namespace=taxon_namespace,
        store_tree_weights=True,
        preserve_underscores=True,
        rooting="default-unrooted",
        ignore_unrecognized_keyword_arguments=True,
    )
    lengths = []
    root_ages = []
    for tree in tree_yielder:
        length = 0.0
        for edge in tree:
            length += edge.length
        lengths.append(length)
        tree.calc_node_ages()
        root_ages.append(tree.seed_node.age)


# not used
# Read data from a CSV file into a NJ_tree
def build_in_nj(filename):  # "distance_matrix.csv"
    '''build in
    nj_tree(is_weighted_edge_distances=True, tree_factory=None)[source]
    Returns an Neighbor-Joining (NJ) tree based on the distances in the matrix.
    Calculates and returns a tree under the Neighbor-Joining algorithm of Saitou and Nei (1987) for the data in the matrix.
    '''
    # Read data from a CSV file into a PhylogeneticDistanceMatrix
    with open(filename) as src:
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            is_allow_new_taxa=True,
            delimiter=",",
        )

    # Calculate the tree
    nj_tree = pdm.nj_tree()
    print(nj_tree.as_string("nexus"))
    return nj_tree


# not used
def build_in_upgma(filename):
    upgma_tree(is_weighted_edge_distances=True, tree_factory=None)[source]
    # Returns an Unweighted Pair Group Method with Arithmetic Mean (UPGMA) tree based on the distances in the matrix.
    # Read data from a CSV file into a PhylogeneticDistanceMatrix
    # object
    with open(filename) as src:
        pdm = dendropy.PhylogeneticDistanceMatrix.from_csv(
            src,
            is_first_row_column_names=True,
            is_first_column_row_names=True,
            is_allow_new_taxa=True,
            delimiter=",",
        )

    # Calculate the tree
    upgma_tree = pdm.upgma_tree()
    print(upgma_tree.as_string("nexus"))
    return upgma_tree


# not used
def read_Distancefile2dict(distance_fpn, name_fp):
    name = open(name_fp, 'r').readlines()
    distance_lst = open(distance_fp, 'r').readlines()
    namelist = []
    for i in range(0, len(name)):
        namelist[i] = float(name[i].strip().split())

    distance_dct = dict((x, dict()) for x in namelist)
    for i in namelist:
        for j in namelist:
            # if i==j:
            #    continue
            distance_dct[i][j] = float(distance_lst[i].strip().split()[j])

    print(distance_dct)
    return distance_dct


# not used
def write_outputfile(fp, tree):
    with open(output_fp, 'w') as f:
        f.write(tree.as_string("newick") + '\n')
    return True


# not used
def write_tree(output_fp, treeparentnodes, treebranchlens):
    treestr_newick = ''
    return treestr_newick


def create_tree(distance_dict, leaves, rootat, alpha):
    tpn, tb = distancedict2treedict_nj(distance_dict)
    print(tpn)
    print(tb)
    tree = TPNandTBL2Tree(tpn, tb, leaves)
    eglst = [eg for eg in tree.postorder_edge_iter()]  # can change the iter method to change order of edges
    rootedge = eglst[rootat]
    rel = rootedge.length
    tree.reroot_at_edge(rootedge, length1=rel * alpha, length2=rel * (1 - alpha))
    print(tree.as_string("newick"))
    print(tree.as_ascii_plot())
    return tree


if "__main__" == __name__:
    test(4, 0.3)

    '''
    parser = OptionParser()
    parser.add_option("-d", "--distance", dest="distance_fp",
                      help="path to the distance file", metavar="FILE")

    parser.add_option("-n", "--name", dest="name_fp",
                      help="path to the species name file", metavar="FILE")

    parser.add_option("-o", "--output", dest = "output_fp",
                      help="path to the output file", metavar="FILE")

    (options, args) = parser.parse_args()
    output_fp = options.output_fp
    distance_dct = read_Distancefile2dict(options.distance_fp,options.name_fp)
    assert(distance_dct!=None)
    tpn,tb=distancedict2treedict_nj(distance_dct)

    '''
