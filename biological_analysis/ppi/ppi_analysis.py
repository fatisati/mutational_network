import pickle as pkl
import time
import networkx as nx
import random
import pandas as pd
import shared_names


def generate_random_network(node_set, node_cnt, link_cnt):
    nodes = random.sample(node_set, node_cnt)
    g = nx.gnm_random_graph(node_cnt, link_cnt)
    mapping = {i: nodes[i] for i in range(len(nodes))}
    g = nx.relabel_nodes(g, mapping)
    return g


class PpiAnalysis:

    def __init__(self, confidence_level, result_path):
        self.res_path = result_path
        print('loading ppi network...')
        self.ppi_network = pkl.load(open(self.res_path + f'ppi-net-{confidence_level}.pkl', 'rb'))
        ppi_net_links = self.ppi_network.edges
        print(f'number of all ppi links in the loaded ppi-network: {len(ppi_net_links)}')

        self.gene_name_dict = pkl.load(open(result_path + shared_names.gene_name_dic, 'rb'))
        # self.vis = vis_utils.Utils(res_path)
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
        for u, v in network.edges:
            g1, g2 = self.gene_name_dict[u], self.gene_name_dict[v]
            if self.has_ppi_interaction(g1, g2):
                ppi_links.append((g1, g2))
        return ppi_links

    def evaluate_coms_ppi_links(self, network: nx.Graph, coms, random_cnt):
        res = []
        all_genes = list(network.nodes)
        enriched_coms_cnt = 0
        min_enrich_score = 1.5
        for com in coms:
            com_net = network.subgraph(com)
            ppi_links = self.get_ppi_links(com_net)
            # print(f'generating {random_cnt} random graphs...')
            st = time.time()
            rand_ppi_cnt = 0
            for _ in range(random_cnt):
                rand_net = generate_random_network(all_genes, len(com), len(com_net.edges))
                rand_ppi_links = self.get_ppi_links(rand_net)
                rand_ppi_cnt += len(rand_ppi_links)
            # print(f'done. took {st - time.time()}')

            rand_ppi_cnt = rand_ppi_cnt / random_cnt
            if rand_ppi_cnt == 0:
                ppi_enrichment = 0
            else:
                ppi_enrichment = len(ppi_links) / rand_ppi_cnt
            res.append({'node-cnt': len(com), 'edge-cnt': len(com_net.edges),
                        'ppi-links': ppi_links, 'ppi-link-cnt': len(ppi_links),
                        'rand-ppi-link_cnt': rand_ppi_cnt,
                        'ppi-enrichment': ppi_enrichment})
            if ppi_enrichment >= min_enrich_score:
                enriched_coms_cnt += 1

        res = pd.DataFrame(res)
        print(f'number of coms with enrichment value equal/greater than {min_enrich_score}: {enriched_coms_cnt}'
              f', ratio to all coms: {enriched_coms_cnt / len(coms)}')
        print(f'average enrichment of all coms: {res["ppi-enrichment"].mean()}')

        all_ppi_links_coms = res[res["ppi-link-cnt"] == res["edge-cnt"]]
        print(f'number of communities which all links are ppi links: {len(all_ppi_links_coms)}')
        # plt.plot(res['size'], res['ppi_link_cnt'], label='original')
        # plt.plot(res['size'], res['rand_ppi_link_cnt'], label='random')
        # plt.xlabel('Community size')
        # plt.ylabel('Number of ppi-mutational links')
        #
        # plt.legend()
        # self.vis.save_and_show_plt('ppi_analysis')
        res.to_excel(self.res_path + f'ppi-analysis-rand{random_cnt}.xlsx')
        return res

    def count_network_ppi_links(self, network):
        ppi_links = self.get_ppi_links(network)
        cnt = len(ppi_links)
        print(f'your network has {cnt} number of ppi links, ratio to all links: {cnt / len(network.edges)}')


if __name__ == '__main__':
    res_path = '../../results/'
    network = pkl.load(open(f'{res_path}network0.15.pkl', 'rb'))
    ppi_analyser = PpiAnalysis(confidence_level='medium', result_path='../../results/')
    ppi_analyser.count_network_ppi_links(network)
    print('----')
    coms = pkl.load(open(res_path + 'h-ensemble1000-coms0.15.pkl', 'rb'))
    print(len(coms))
    ppi_analyser.evaluate_coms_ppi_links(network, coms, 1000)
