import pandas as pd

import shared_names


class GeneLoc:
    def __init__(self, result_folder):
        self.result_folder = result_folder
        self.gene_loc_dic = pkl.load(open(result_folder + shared_names.gene_location_dict, 'rb'))

    def distance(self, g1, g2):
        if g1 not in self.gene_loc_dic:
            print(f'{g1} not found in gene-loc')
            return -1
        if g2 not in self.gene_loc_dic:
            print(f'{g2} not found in gene-loc')
            return -1
        loc1 = self.gene_loc_dic[g1]['mid']
        loc2 = self.gene_loc_dic[g2]['mid']

        return abs(loc2 - loc1)

    def check_dic_contain(self, gene):
        if gene not in self.gene_loc_dic:
            print(f'{gene} not found in gene-loc dic')
            return False
        return True

    def analyse_cis_trans_distance_strand(self, com_net):
        cis = []
        trans = []
        distance_sum = 0
        same_strand = []

        for g1, g2 in com_net.edges:
            if not (self.check_dic_contain(g1)):
                continue
            if not (self.check_dic_contain(g2)):
                continue

            if self.gene_loc_dic[g1]['chr'] == self.gene_loc_dic[g2]['chr']:
                cis.append((g1, g2))
                distance_sum += self.distance(g1, g2)
            else:
                trans.append((g1, g2))

            if self.gene_loc_dic[g1]['strand'] == self.gene_loc_dic[g2]['strand']:
                same_strand.append((g1, g2))
        if len(cis) > 0:
            avg_dist = distance_sum / len(cis)
        else:
            avg_dist = 0
        return cis, trans, avg_dist, same_strand

    def analyze_all_coms(self, network, coms):
        cis, trans, dist, strand = self.analyse_cis_trans_distance_strand(network)
        link_cnt = len(network.edges)
        print(
            f'the network has {len(cis) / link_cnt} cis and {len(trans) / link_cnt} trans links (ratio to all links). '
            f'the avg distance between genes of cis links: {dist}, ratio of links with genes on the same strand: {len(strand) / link_cnt}')

        subnet_list = [network.subgraph(com) for com in coms]
        cis_list, trans_list, distance_list, strand_list = [], [], [], []

        for subnet in subnet_list:
            cis, trans, dist, strand = self.analyse_cis_trans_distance_strand(subnet)
            cis_list.append(cis)
            trans_list.append(trans)
            distance_list.append(dist)
            strand_list.append(strand)
        save_path = self.result_folder + 'gene-location-analysis.csv'
        pd.DataFrame({'cis': cis_list, 'trans': trans_list, 'avg-cis-distance': distance_list,
                      'same-strand': strand_list}).to_csv(save_path)
        return cis_list, trans_list, distance_list, strand_list


import pickle as pkl

if __name__ == '__main__':
    res_path = '../results/'
    coms = pkl.load(open(res_path + 'h-ensemble1000-coms0.15.pkl', 'rb'))
    network = pkl.load(open(res_path + 'network0.15.pkl', 'rb'))
    gene_loc_analyzer = GeneLoc(res_path)
    gene_loc_analyzer.analyze_all_coms(network, coms)
