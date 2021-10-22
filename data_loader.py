import pickle as pkl

import networkx as nx
import pandas as pd
import config

import config
from utils import *
from biological_analysis.ppi_interactions import *


def load_name_coms(thr):
    coms = load_coms(thr)
    translate_com = lambda com: [Dictionary().get_gene_name(gene) for gene in com]
    return [translate_com(com) for com in coms]


def load_coms(thr):
    return pkl.load(open(config.community_path + 'coms' + str(thr) + '.pkl', 'rb'))


def load_name_network(thr):
    net = load_network(thr)
    nodes = list(net.nodes)
    new_net = nx.relabel_nodes(net, lambda x: Dictionary().get_gene_name(x))
    return new_net


def load_network(thr):
    print(f'loading network thr: {thr}')
    net = pkl.load(open(config.network_path + f'network{thr}.pkl', 'rb'))
    print('---done---')
    return net


def load_gene_name_dic():
    print('loading gene-name dic')
    dic = pkl.load(open(config.gene_name_dic_path, 'rb'))
    print('---done---')
    return dic


def load_ppi_interaction_network(min_score):
    print('loading interaction network')
    confidence = get_confidence_score(min_score)
    file_path = f'{config.stringdb_path}interaction_network_{confidence}.pkl'

    try:
        network = pkl.load(open(file_path, 'rb'))
    except:
        print('network not found. generating network...')
        interactions = load_ppi_interaction_dic(min_score)
        network = nx.Graph(interactions)
        pkl.dump(network, open(file_path, 'wb'))

    prots = set(network.nodes)
    name = Dictionary().map_prot_name(list(prots)[0])
    print(
        f'number of proteins: {len(prots)}, number of proteins with available mappings: {len(prots.intersection(set(Dictionary().prot_gene.keys())))}')
    network = nx.relabel_nodes(network, Dictionary().map_prot_name)

    return network


def load_ppi_interaction_dic(min_score):
    confidence = get_confidence_score(min_score)
    path = f'{config.dump_path}stringdb_ppi_{confidence}.pkl'

    try:
        interactions_dic = pkl.load(open(path, 'rb'))
    except:
        interactions_dic = generate_ppi_dic(min_score)

    return interactions_dic


def load_gene_prot_dics():
    gene_prot_path = config.dump_path + 'gene_prot.pkl'
    prot_gene_path = config.dump_path + 'prot_gene.pkl'
    # gene_name_path = util.res_path()+'gene_name.pkl'

    gene_prot = pkl.load(open(gene_prot_path, 'rb'))
    prot_gene = pkl.load(open(prot_gene_path, 'rb'))
    # gene_name = pkl.load(open(gene_name_path, 'rb'))

    return gene_prot, prot_gene


def load_table(name):
    return pd.read_excel(f'{config.table_path}{name}').fillna('')


def load_coms_ontology():
    return load_table(config.coms_enrichment_name)[['Process', 'Component', 'Function']]


def load_coms_pathways():
    return load_table(config.coms_pathways)


def load_cancer_pathways_genes():
    return pd.read_excel(config.cancer_pathways_genes_name)


def load_pmids():
    return load_table('cancer_genes_pmid.xlsx'), load_table('breast_cancer_genes_pmid.xlsx')


def load_breast_cancer_genes():
    cancer_gene_df = load_table('breast_cancer_genes_pmid')
    cancer_genes = set(cancer_gene_df[cancer_gene_df['breast cancer'] != '']['gene'])
    return cancer_genes


def load_cancer_genes():
    cancer_gene_df = load_table('cancer_genes_pmid')
    cancer_genes = set(cancer_gene_df[cancer_gene_df['cancer pmid'] != '']['gene'])
    print(len(cancer_genes))
    return cancer_genes


def load_gene_chromosome():
    df = pd.read_excel(config.data_path + 'genes_list_hg19.xlsx')
    return df


def load_gene_patient():
    path = config.old_dump_path + 'gene_patient.pkl'
    return pkl.load(open(path, 'rb'))


def load_patient_gene():
    path = config.old_dump_path + 'patient_gene.pkl'
    return pkl.load(open(path, 'rb'))


def load_assortativity_table():

    return pd.read_excel(config.old_table_path + 'assortativity.xlsx', index_col=0, header=None)
