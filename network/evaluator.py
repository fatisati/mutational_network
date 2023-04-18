import pandas as pd
import pickle as pkl


class NetworkEvaluator:
    def __init__(self, network):
        self.network = network
        self.gene_table = pd.read_excel('../../data/genes_list_hg19.xlsx')
        coding_types = ['protein_coding', 'processed_transcript']
        coding_rows = self.gene_table[self.gene_table['Gene_class'].isin(coding_types)]
        self.coding_genes = list(coding_rows['gene_name'])

    def evaluate(self):
        node_cnt = len(self.network.nodes)
        edge_cnt = len(self.network.edges)

        coding_nodes = set(self.coding_genes).intersection(self.network.nodes)
        non_coding_nodes = set(self.network.nodes) - set(coding_nodes)

        coding_links = self.network.subgraph(coding_nodes).edges
        non_coding_links = self.network.subgraph(non_coding_nodes).edges

        print(f'number of nodes: {node_cnt}., number of edges: {edge_cnt}')
        print(f'number of coding genes: {len(coding_nodes)}, non-coding: {len(non_coding_nodes)}')
        print(f'number of cc links:: {len(coding_links)}, nn links: {len(non_coding_links)}')


if __name__ == '__main__':
    network = pkl.load(open('../results/network.15.pkl', 'rb'))
    evaluator = NetworkEvaluator(network)
    evaluator.evaluate()
