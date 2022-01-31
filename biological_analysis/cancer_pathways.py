class CancerPathways:
    def __init__(self, data_loader):
        self.cancer_pathway_genes = data_loader.load_cancer_pathways_genes()

    def find_com_cancer_pathways(self, com):
        com_pathways = []
        for _, row in self.cancer_pathway_genes.iterrows():
            common_genes = set(com).intersection(set(row['genes']))
            if len(common_genes) > 0:
                com_pathways.append(row['pathways'])
        return com_pathways


if __name__ == '__main__':
    cp = CancerPathways('../../../data/')
    sample = cp.cancer_pathway_genes['genes'].iloc[0]
    print(sample)
    print(type(sample))
