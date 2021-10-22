import pandas as pd

import data_loader


class Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]


class Dictionary(metaclass=Singleton):
    def __init__(self):
        self.gene_name = data_loader.load_gene_name_dic()
        self.gene_prot, self.prot_gene = None, None
        self.gene_chromosome = None

        print('---gene name dictionary loaded---')

    def get_gene_name(self, gene_id):
        if gene_id in self.gene_name:
            return self.gene_name[gene_id]['name']
        return gene_id

    def get_coms_names(self, coms):
        coms_names = []
        for com in coms:
            coms_names.append(map(self.get_gene_name, com))
        return coms_names

    def get_prot_gene(self, prot):

        if not self.prot_gene:
            print('loading prot-gene dic')
            self.gene_prot, self.prot_gene = data_loader.load_gene_prot_dics()

        return get_if_has(self.prot_gene, prot)

    def map_prot_name(self, prot):

        gene_id = self.get_prot_gene(prot)

        gene_name = self.get_gene_name(gene_id)
        return gene_name

    def get_gene_chromosome(self, gene_name):
        if type(self.gene_chromosome) != type(pd.DataFrame([])):
            print(f'loading gene chrome...')
            self.gene_chromosome = data_loader.load_gene_chromosome()
            print('---done---')

        try:
            return self.gene_chromosome[self.gene_chromosome['gene_symbol'] == f'"{gene_name}"']['chr'].iloc[0]
        except:
            print(f'gene-chrome not found for gene {gene_name}')


def add_dic_set(dic, k, v, score):
    if k not in dic:
        dic[k] = {}
    dic[k][v] = score


def get_if_has(dic, k):
    try:
        return dic[k]
    except:
        return 'not-found-' + k


def find_com(coms, gene):
    for com in coms:
        if gene in com:
            return com
    return 'not found'
