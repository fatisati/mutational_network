from networkx.algorithms import community
import community

from ecg import *

from data_loader import *
from save_data import save_coms


def community_list_dic(com_dic):
    ans_dic = {}

    for gin in com_dic:
        if com_dic[gin] not in ans_dic:
            ans_dic[com_dic[gin]] = []
        ans_dic[com_dic[gin]].append(gin)

    ans = []
    for k in ans_dic.keys():
        ans.append(ans_dic[k])

    # print(ans)
    return ans


def louvain(g):
    return community_list_dic(community.best_partition(g))


def consensus(g):
    g2 = ig.Graph.TupleList(g.edges())
    coms = g2.community_ecg()

    coms_list = []
    for item in coms:
        com = []
        for gene in item:
            com.append(g2.vs[gene]['name'])
        coms_list.append(com)
    return coms_list


def hierarchical(network, community_detection_func):
    coms = community_detection_func(network)
    if len(coms) < 2:
        return coms

    res = []
    for com in coms:
        subgraph = network.subgraph(com)
        deep_coms = hierarchical(subgraph, community_detection_func)
        for dcom in deep_coms:
            res.append(dcom)
    return res


def hierarchical_ensemble(g):
    print('ensemble hierarchical')
    coms = hierarchical(g, consensus)
    print('---done---')
    return coms


def hierarchical_louvain(g):
    print('louvain hierarchical')
    coms = hierarchical(g, louvain)
    print('---done---')
    return coms


if __name__ == '__main__':
    selected_threshold = 0.15
    experiment_count = 3

    network = load_network(selected_threshold)

    algos = {'louvain': louvain, 'ensemble': consensus,
             'h_louvain': hierarchical_louvain, 'h_ensemble': hierarchical_ensemble}
    for i in range(experiment_count):
        for alg_name in algos:
            coms = algos[alg_name](network)
            save_coms(coms, f'{alg_name}_run{i + 3}')
