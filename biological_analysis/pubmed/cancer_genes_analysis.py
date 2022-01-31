import pandas as pd


class CancerGene:
    def __init__(self, data_folder):
        self.cancer_genes = pd.read_excel(data_folder + 'pubmed/cancer_genes.xlsx')
        self.breast_cancer_genes = pd.read_excel(data_folder + 'pubmed/breast_cancer_genes.xlsx')

    def com_cancer_breast_cancer_genes(self, com):
        com_cancer = set(self.cancer_genes['gene']).intersection(set(com))
        com_breast_cancer = set(self.breast_cancer_genes['gene']).intersection(set(com))
        return com_cancer, com_breast_cancer
