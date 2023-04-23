from networkx.algorithms import community
import community

from community_.ecg import *
import networkx as nx
import pickle as pkl


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


def consensus(g, ens_size=100):
    edges, weights = zip(*nx.get_edge_attributes(g, 'weight').items())
    g2 = ig.Graph.TupleList(g.edges())
    # for e in g2.es:
    #     print(e.source, e.target)
    coms = g2.community_ecg(weights=weights, ens_size=ens_size)

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
        # deep_coms_sizes = [len(c) for c in deep_coms]
        # print(f'{len(com)} -> {deep_coms_sizes}')
        for dcom in deep_coms:
            res.append(dcom)
    return res


def hierarchical_ensemble(g, ens_size=100):
    print('ensemble hierarchical')
    ensemble_func = lambda g: consensus(g, ens_size)
    coms = hierarchical(g, ensemble_func)
    print('---done---')
    return coms


def hierarchical_louvain(g):
    print('louvain hierarchical')
    coms = hierarchical(g, louvain)
    print('---done---')
    return coms


if __name__ == '__main__':
    result_path = '../results/'
    network = pkl.load(open(f'{result_path}/network0.15.pkl', 'rb'))
    print(len(network.nodes), len(network.edges))

    ensemble_size = 1000
    coms = hierarchical_ensemble(network, ensemble_size)
    print(f'number of founded communities: {len(coms)}')
    save_name = f'h-ensemble{ensemble_size}-coms0.15'
    pkl.dump(coms, open(f'{result_path}/{save_name}.pkl', 'wb'))
