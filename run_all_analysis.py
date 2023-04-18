import pandas as pd
from data_loader import DataLoader
from community import community_detection
from community_analysis import CommunityAnalysis
from biological_analysis.ppi.ppi_analysis import PpiAnalysis
from biological_analysis.string_api import StringApi
from biological_analysis.cancer_pathways import CancerPathways
from biological_analysis.gene_locations import GeneLoc
from biological_analysis.pubmed.cancer_genes_analysis import CancerGene


def set_col_arr_and_arr_cnt(df, col, col_values):
    df[col] = col_values
    df[f'{col}_cnt'] = [len(arr) for arr in col_values]
    return df


if __name__ == '__main__':
    data_folder = '../../data/'
    save_path = '../../results/new-run/'
    data_loader = DataLoader(data_folder, save_path)

    # network_thr = 0.15
    # network = data_loader.load_network(network_thr)

    # all_analysis_df = pd.DataFrame([])

    # ensemble_size = 5
    # coms = community.hierarchical_ensemble(network, ensemble_size)
    # data_loader.save_com(0.15, coms)
    # all_analysis_df['coms'] = coms
    coms = data_loader.load_com(0.15)

    # com_analysis = CommunityAnalysis(coms, save_path, data_loader)
    # coms_coverages = com_analysis.sample_coverages()
    # all_analysis_df['patient_coverages'] = coms_coverages
    # all_analysis_df.to_excel(save_path + 'all_analysis.xlsx')

    # rand_cnt = 10
    # ppi_analysis = PpiAnalysis(data_folder, save_path)
    # ppi_df = ppi_analysis.evaluate_coms_ppi_links(network, coms, rand_cnt)
    # all_analysis_df = pd.concat([all_analysis_df, ppi_df], axis = 1)
    # all_analysis_df.to_excel(save_path + 'all_analysis.xlsx')

    # string_api = StringApi()
    # all_kegg, all_go = [], []
    # for com in coms:
    #     kegg, go = string_api.get_enrichment(com)
    #     all_kegg.append(kegg)
    #     all_go.append(go)
    # all_analysis_df = set_col_arr_and_arr_cnt(all_analysis_df, 'kegg', all_kegg)
    # all_analysis_df = set_col_arr_and_arr_cnt(all_analysis_df, 'go', all_go)
    # all_analysis_df.to_excel(save_path + 'all_analysis.xlsx')

    cancer_pathways = CancerPathways(data_folder)
    coms_cancer_pathways = [cancer_pathways.find_com_cancer_pathways(com) for com in coms]
    # all_analysis_df = set_col_arr_and_arr_cnt(all_analysis_df, 'cancer_pathways', coms_cancer_pathways)
    # all_analysis_df.to_excel(save_path + 'all_analysis.xlsx')

    # gene_loc = GeneLoc(data_folder)
    # cis_trans_dist_strand = {'cis': [], 'trans': [], 'dist': [], 'strand': []}
    # for com in coms:
    #     com_net = network.subgraph(com)
    #     res = gene_loc.analyse_cis_trans_distance_strand(com_net)
    #     cnt = 0
    #     for key in cis_trans_dist_strand.keys():
    #         cis_trans_dist_strand[key].append(res[cnt])
    #         cnt += 1
    # all_analysis_df['avg_cis_distance'] = cis_trans_dist_strand['dist']
    # for key in ['cis', 'trans', 'strand']:
    #     all_analysis_df[f'{key}_links'] = cis_trans_dist_strand[key]
    #     all_analysis_df[f'{key}_cnt'] = [len(arr) for arr in cis_trans_dist_strand[key]]
    # all_analysis_df.to_excel(save_path + 'all_analysis.xlsx')
    #
    # cancer_gene = CancerGene(data_folder)
    # coms_cancer_genes, coms_breast_cancer_genes = [], []
    # for com in coms:
    #     caner, breast = cancer_gene.com_cancer_breast_cancer_genes(com)
    #     coms_cancer_genes.append(caner)
    #     coms_breast_cancer_genes.append(breast)
    # all_analysis_df = set_col_arr_and_arr_cnt(all_analysis_df, 'cancer_genes', coms_cancer_genes)
    # all_analysis_df = set_col_arr_and_arr_cnt(all_analysis_df, 'breast_cancer_genes', coms_breast_cancer_genes)
    #
    # all_analysis_df.to_excel(save_path + 'all_analysis.xlsx')
