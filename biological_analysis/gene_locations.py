import data_loader
import utils
from biological_analysis.network_analysis import *
from utils import *
from copy import deepcopy
import pandas as pd
import config


def get_gene_chrome_dic(gene, com, network):
    neighbors = find_neighbors(gene, network)
    print(f'neighbor count: {len(neighbors)}')

    chrome_dic = {'same chrome': [], 'diff chrome': []}
    res_dic = {'neighbor': deepcopy(chrome_dic),
               'same com': deepcopy(chrome_dic),
               'same com and neighbor': deepcopy(chrome_dic)}

    dict = Dictionary()
    gene_chrome = dict.get_gene_chromosome(gene)
    print(f'gene chrome: {gene_chrome}')
    res_dic['chromosome'] = gene_chrome
    def add_gene_dic(check_gene, gene_relation):
        if check_gene == gene:
            return
        if dict.get_gene_chromosome(check_gene) == gene_chrome:
            chrome_label = 'same chrome'
        else:
            chrome_label = 'diff chrome'
        res_dic[gene_relation][chrome_label].append(check_gene)

    for neighbor in neighbors:
        add_gene_dic(neighbor, 'neighbor')

    print(f'number of genes in same com: {len(com)}')
    for com_gene in com:
        add_gene_dic(com_gene, 'same com')

    com_links = set(com).intersection(set(neighbors))

    for com_link in com_links:
        add_gene_dic(com_link, 'same com and neighbor')
    return res_dic


def analyze_genes_chrome(genes, coms, network):
    res_df = []
    for gene in genes:
        com = utils.find_com(coms, gene)

        chr_dic = get_gene_chrome_dic(gene, com, network)
        chr_dic['gene'] = gene
        res_df.append(chr_dic)
    pd.DataFrame(res_df).to_excel(config.save_path + 'gene-chromes.xlsx')


if __name__ == '__main__':
    coms = data_loader.load_name_coms(0.15)
    network = data_loader.load_name_network(0.15)

    genes = ['MAGI2']
    analyze_genes_chrome(genes, coms, network)
