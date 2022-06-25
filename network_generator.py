from scipy.sparse import csr_matrix
import numpy as np
import networkx as nx


class MutationalNetwork:
    def __init__(self, data_folder):
        self.data_folder = data_folder
        self.data_path = data_folder + 'sample.tsv'
        self.gene_key = 'gene_affected'
        self.patient_key = 'icgc_specimen_id'
        self.gene_patient = {}
        self.all_genes = set([])
        self.all_patient = set([])

    def generate_gene_patient(self):
        data = open(self.data_path)
        header = data.readline().split('\t')
        gene_idx = header.index(self.gene_key)
        patient_idx = header.index(self.patient_key)

        row = data.readline()
        while row:
            row = row.split('\t')
            self.add_patient_to_gene(row[gene_idx], row[patient_idx])
            row = data.readline()

        print(f'gene patient created. patient count: {len(self.all_patient)}, gene cnt: {len(self.all_genes)}')

    def add_patient_to_gene(self, gene, patient):
        if gene not in self.gene_patient:
            self.gene_patient[gene] = set([])
            self.all_genes.add(gene)
        self.gene_patient[gene].add(patient)
        self.all_patient.add(patient)

    def gene_patient_mat(self):
        mat = []
        for gene in self.all_genes:
            gene_row = []
            for patient in self.all_patient:
                exist = (patient in self.gene_patient[gene])
                gene_row.append(int(exist))
            mat.append(gene_row)
        return np.array(mat)

    def generate_network(self, threshold):
        self.generate_gene_patient()
        mat = self.gene_patient_mat()
        transposed = mat.T

        gene_gene = (csr_matrix(mat) * csr_matrix(transposed)).toarray()
        print(f'gene-gene: {gene_gene}, min: {gene_gene.min()}, max: {gene_gene.max()}')

        graph = nx.Graph()
        genes = list(self.all_genes)
        for i in range(len(self.all_genes)):
            for j in range(i, len(self.all_genes)):
                if gene_gene[i][j] >= threshold * len(self.all_patient):
                    graph.add_edge(genes[i], genes[j], weight=gene_gene[i][j])
        print(f'graph edge cnt: {len(graph.edges)}, node cnt: {len(graph.nodes)}')
        return graph

    def save_network(self):
        pass


def make_sample_data(data_folder, n_rows=1000):
    data_path = data_folder + 'simple_somatic_mutation.open.tsv'
    data = open(data_path)
    sample = open(data_folder + 'sample.tsv', 'w')

    header = data.readline()
    sample.write(header)

    for _ in range(n_rows):
        line = data.readline()
        sample.write(line)
    sample.close()


if __name__ == '__main__':
    # make_sample_data('../../data/', 1000000)
    MutationalNetwork('../../data/').generate_network(0.1)
