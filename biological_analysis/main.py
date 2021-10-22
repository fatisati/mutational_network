import data_loader
import numpy as np

from ppi_interactions import *
from gene_locations import *

def select_most_related_brest_cancer_genes(n):
    cancer, breast_cancer = data_loader.load_pmids()
    pmid_cnts = []
    for index, row in breast_cancer.iterrows():
        pmid_cnts.append(len(row['breast cancer'].split(',')))

    selected_genes_idx = np.array(pmid_cnts).argsort()[-n:]
    selected_genes_idx = selected_genes_idx[::-1]

    selected_genes = []
    for idx in selected_genes_idx:
        print(breast_cancer.iloc[idx]['gene'], pmid_cnts[idx])
        selected_genes.append(breast_cancer.iloc[idx]['gene'])
    return selected_genes

if __name__ == '__main__':
    selected_genes = select_most_related_brest_cancer_genes(10)
    network = data_loader.load_name_network(0.15)
    coms = data_loader.load_name_coms(0.15)

    analyze_genes_chrome(selected_genes, coms, network)
    # ppi = PpiAnalysis()
    # for gene in selected_genes:
    #     ppi.analyze_gene(gene, network)