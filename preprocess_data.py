import pickle as pkl
import shared_names
import pandas as pd


class DataPreprocessor:
    def __init__(self, data_folder, result_folder, mutation_data='sample-mutations.tsv'):
        self.data_folder = data_folder
        self.mutation_path = data_folder + mutation_data
        self.result_folder = result_folder
        self.gene_key = 'gene_affected'
        self.patient_key = 'icgc_specimen_id'
        self.gene_patient = {}
        self.all_genes = set([])
        self.all_patient = set([])

    def generate_gene_patient(self):
        data = open(self.mutation_path)
        header = data.readline().split('\t')
        gene_idx = header.index(self.gene_key)
        patient_idx = header.index(self.patient_key)

        row = data.readline()
        while row:
            row = row.split('\t')
            self.add_patient_to_gene(row[gene_idx], row[patient_idx])
            row = data.readline()
        print(f'gene-patient created. patient count: {len(self.all_patient)}, gene cnt: {len(self.all_genes)}')
        pkl.dump(self.gene_patient, open(f'{self.result_folder}/{shared_names.gene_patient}', 'wb'))

    def generate_patient_gene(self):
        if len(self.gene_patient) == 0:
            self.generate_gene_patient()
        patient_gene = {}
        for gene in self.gene_patient:
            for patient in self.gene_patient[gene]:
                if patient not in patient_gene:
                    patient_gene[patient] = []
                patient_gene[patient].append(gene)
        pkl.dump(patient_gene, open(f'{self.result_folder}/{shared_names.patient_gene}', 'wb'))

    def add_patient_to_gene(self, gene, patient):
        if len(gene) == 0:
            return
        if gene not in self.gene_patient:
            self.gene_patient[gene] = set([])
            self.all_genes.add(gene)
        self.gene_patient[gene].add(patient)
        self.all_patient.add(patient)

    def generate_gene_name_dict(self):
        gene_list = pd.read_excel(self.data_folder + 'genes_list_hg19.xlsx')

        gene_name_dict = {}
        for index, row in gene_list.iterrows():
            gene_id = row['gene_name']
            gene_name = row['gene_symbol']
            if gene_id in gene_name_dict:
                print(
                    f'duplicate names for gene-id: {gene_id}. new-name: {gene_name}, prev-name: {gene_name_dict[gene_id]}')
                continue
            # remove double-quote ""
            gene_name_dict[gene_id] = gene_name[1:-1]
        save_path = self.result_folder + shared_names.gene_name_dic
        pkl.dump(gene_name_dict, open(save_path, 'wb'))
        print(f'gene-name-dict created and saved in {save_path}. number of keys: {len(gene_name_dict)}')


def make_sample_data(data_folder, n_rows=1000):
    data_path = data_folder + 'simple_somatic_mutation.open.tsv'
    data = open(data_path)
    sample = open(data_folder + 'sample-mutations.tsv', 'w')

    header = data.readline()
    sample.write(header)

    for _ in range(n_rows):
        line = data.readline()
        sample.write(line)
    sample.close()


if __name__ == '__main__':
    data_folder = './data/'
    data_name = 'simple_somatic_mutation.open.tsv'
    result_folder = './results/'

    data_preprocessor = DataPreprocessor(data_folder, result_folder, data_name)
    # data_preprocessor.generate_gene_patient()
    # data_preprocessor.generate_patient_gene()
    data_preprocessor.generate_gene_name_dict()
