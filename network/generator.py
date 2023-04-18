from scipy.sparse import csr_matrix
import numpy as np
import networkx as nx
import pickle as pkl
import shared_names


class MutationalNetwork:
    def __init__(self, result_folder):
        self.data_folder = result_folder
        self.gene_patient_path = result_folder + shared_names.gene_patient
        self.gene_patient, self.all_genes, self.all_patient = self.load_gene_patient()
        self.out_folder = result_folder

    def load_gene_patient(self):
        gene_patient = pkl.load(open(self.gene_patient_path, 'rb'))
        all_genes = list(gene_patient.keys())
        all_patients = set([])
        for gene in all_genes:
            for patient in gene_patient[gene]:
                all_patients.add(patient)
        print(
            f'{self.gene_patient_path} loaded. number of genes: {len(all_genes)}, number of patient: {len(all_patients)}')
        return gene_patient, all_genes, list(all_patients)

    def remove_rare_genes(self, rare_gene_thr):
        patient_cnt = len(self.all_patient)
        filtered_genes = set([])
        for gene in self.all_genes:
            mutation_ratio = len(self.gene_patient[gene]) / patient_cnt
            if mutation_ratio >= rare_gene_thr:
                filtered_genes.add(gene)
        print(f'number of genese after rare removal: {len(filtered_genes)}')
        return filtered_genes

    def gene_patient_mat(self, genes):
        mat = []
        for gene in genes:
            gene_row = []
            for patient in self.all_patient:
                exist = (patient in self.gene_patient[gene])
                gene_row.append(int(exist))
            mat.append(gene_row)
        return np.array(mat)

    def generate_network(self, threshold):
        filtered_genes = self.remove_rare_genes(threshold)
        mat = self.gene_patient_mat(filtered_genes)

        gene_gene = (csr_matrix(mat) * csr_matrix(mat).transpose()).toarray()
        print(f'gene-gene: min: {gene_gene.min()}, max: {gene_gene.max()}')

        graph = nx.Graph()
        genes = list(filtered_genes)
        all_patient_cnt = len(self.all_patient)
        for i in range(len(filtered_genes)):
            for j in range(i + 1, len(filtered_genes)):
                if gene_gene[i][j] >= threshold * all_patient_cnt:
                    graph.add_edge(genes[i], genes[j], weight=gene_gene[i][j] / all_patient_cnt)
        print(f'graph edge cnt: {len(graph.edges)}, node cnt: {len(graph.nodes)}')
        return graph

    def save_network(self, graph, name):
        save_path = f'{self.out_folder}/{name}.pkl'
        pkl.dump(graph, open(save_path, 'wb'))
        print(f'graph saved in {save_path}')


if __name__ == '__main__':
    mutnet = MutationalNetwork('../results/')
    graph = mutnet.generate_network(0.15)
    mutnet.save_network(graph, 'network0.15')
#
