import pickle as pkl
import time

import matplotlib.pyplot as plt
import networkx as nx
import random
import pandas as pd
import visualization


def generate_random_network(node_set, node_cnt, link_cnt):
    nodes = random.sample(node_set, node_cnt)
    g = nx.gnm_random_graph(node_cnt, link_cnt)
    mapping = {i: nodes[i] for i in range(len(nodes))}
    g = nx.relabel_nodes(g, mapping)
    return g


class PpiAnalysis:

    def __init__(self, data_folder, res_path):
        self.string_folder = data_folder + '/stringdb/'
        self.res_path = res_path
        print('loading ppi network...')
        self.ppi_network = pkl.load(open(self.string_folder + 'interaction_network_medium.pkl', 'rb'))
        print(self.ppi_network.nodes)
        print('done')
        self.vis = vis_utils.Utils(res_path)
        self.min_path_length = 2

    def has_ppi_interaction(self, p1, p2):
        try:
            path_len = nx.shortest_path_length(self.ppi_network, source=p1, target=p2)
        except Exception as e:
            return False
        if path_len <= self.min_path_length:
            return True
        return False

    def get_ppi_links(self, network: nx.Graph):

        ppi_links = []
        for p1, p2 in network.edges:
            if self.has_ppi_interaction(p1, p2):
                ppi_links.append((p1, p2))
        return ppi_links

    def evaluate_coms_ppi_links(self, network: nx.Graph, coms, random_cnt):
        res = []
        all_genes = list(network.nodes)
        for com in coms:
            com_net = network.subgraph(com)
            ppi_links = self.get_ppi_links(com_net)
            print(f'number of ppi links: {len(ppi_links)}, edge cnt: {len(com_net.edges)}')

            print(f'generating {random_cnt} random graphs...')
            st = time.time()
            rand_ppi_cnt = 0
            for _ in range(random_cnt):
                rand_net = generate_random_network(all_genes, len(com), len(com_net.edges))
                rand_ppi_links = self.get_ppi_links(rand_net)
                rand_ppi_cnt += len(rand_ppi_links)
            print(f'done. took {st - time.time()}')

            rand_ppi_cnt = rand_ppi_cnt / random_cnt
            if rand_ppi_cnt == 0:
                ppi_enrichment= 0
            else:
                ppi_enrichment = len(ppi_links) / rand_ppi_cnt
            res.append({'com': com, 'size': len(com),
                        'ppi_links': ppi_links, 'ppi_link_cnt': len(ppi_links),
                        'rand_ppi_link_cnt': rand_ppi_cnt,
                        'ppi_enrichment': ppi_enrichment})

        res = pd.DataFrame(res)
        res.to_excel(self.res_path + f'ppi_analysis_rand{random_cnt}.xlsx')

        plt.plot(res['size'], res['ppi_link_cnt'], label='original')
        plt.plot(res['size'], res['rand_ppi_link_cnt'], label='random')
        plt.xlabel('Community size')
        plt.ylabel('Number of ppi-mutational links')

        plt.legend()
        self.vis.save_and_show_plt('ppi_analysis')
        return res
