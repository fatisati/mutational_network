import utils
import matplotlib.pyplot as plt
import seaborn as sns
from data_loader import DataLoader

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

    def show(self, coms):
        coms_pathways = [self.find_com_cancer_pathways(com) for com in coms]
        x, y, data = utils.make_heatmap_data(coms_pathways)
        sns.heatmap(data)
        plt.show()

if __name__ == '__main__':
    data_loader = DataLoader('/home/fatemeh/my-projects/mutational-network/', '')
    cp = CancerPathways(data_loader)
    coms = data_loader.load_com(0.15)
    cp.show(coms)
