import config
import data_loader
import utils
from biological_analysis.network_analysis import *
import pickle as pkl


def get_confidence_score(min_score):
    if min_score <= 400:
        return 'low'
    elif min_score <= 700:
        return 'medium'

    return 'high'


def generate_ppi_dic(min_score=0):
    pass


class PpiAnalysis:
    def __init__(self):
        self.ppi_network = data_loader.load_ppi_interaction_network(401)

    def analyze_gene(self, gene, network):
        print(gene)
        gene_ppi_interactions = list(self.ppi_network[gene])
        neighbors = find_neighbors(gene, network)

        ppi_mutational = set(neighbors).intersection(set(gene_ppi_interactions))

        ppi_mutational_ratio = len(ppi_mutational) / len(gene_ppi_interactions)
        mutational_ppi_ratio = len(ppi_mutational) / len(neighbors)

        print(
            f'this gene has {len(gene_ppi_interactions)} ppi interactions.')

        if round(ppi_mutational_ratio, 2) > 0:
            print(f'{int(ppi_mutational_ratio*100)} percent of them are in mutational network too')

        print(
            f'this gene has {len(neighbors)} mutational interactions.')

        if round(mutational_ppi_ratio, 2) > 0:
            print(f'{int(mutational_ppi_ratio*100)} percent of them are ppi interaction too.')

        if len(ppi_mutational) > 0:
            print(f'number of mutational and ppi interaction links {len(ppi_mutational)}')
            print('ppi-mutational links:')
            for neighbor in ppi_mutational:
                weight = network[gene][neighbor]["weight"]
                print(f'{neighbor}: co-mutation rate ({round(weight, 2)})')

        # com_neighbors = set(com).intersection(set(neighbors))
        # print(
        #     f'this gene has {len(com_neighbors)} neighbors in its community which is {len(com_neighbors) / len(neighbors)} percent of all neigbors')
        # com_ppi = set(gene_ppi_interactions).intersection(com)
        # print(
        #     f'this gene has {len(com_ppi)} ppi interactions in its community which is {len(com_ppi) / len(gene_ppi_interactions)} percent of all gene ppi interactions')
        #
        # com_ppi_mutational = set(ppi_mutational).intersection(com)
        # print(f'this gene has {len(com_ppi_mutational)} ppi-mutational in its community')
        # print(com_ppi_mutational)
        # print('-------------------')

if __name__ == '__main__':
    gene = 'MAGI2'
    thr = 0.15
    network = data_loader.load_name_network(thr)
    coms = data_loader.load_coms(0.15)

    com = find_com(coms, gene)
    ppi = PpiAnalysis()
    ppi.analyze_gene(gene, network, com)
