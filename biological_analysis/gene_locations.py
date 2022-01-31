import pandas as pd


class GeneLoc:
    def __init__(self, data_folder):
        loc_df = pd.read_excel(data_folder + 'genes_list_hg19.xlsx')
        self.gene_loc_dic = {}
        print('generating gene-loc dic...')
        for _, row in loc_df.iterrows():
            gene = row['gene_symbol'][1:-1]
            chr = row['chr']
            start = int(row['start'])
            end = int(row['end'])
            strand = row['strand']
            gene_id = row['gene_name']

            dic_row = {'chr': chr, 'start': start, 'end': end,
                       'mid': (start + end) / 2, 'strand': strand}

            self.gene_loc_dic[gene] = dic_row
            self.gene_loc_dic[gene_id] = dic_row
        print('done')

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

        return cis, trans, distance_sum / len(cis), same_strand
