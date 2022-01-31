import pandas as pd
import pickle as pkl


class GeneNameDic:
    def __init__(self, data_folder):
        self.data_folder = data_folder
        self.dic = pkl.load(open(data_folder + 'dump/gene_name.pkl', 'rb'))

    def generate_dic(self):
        print('generating gene-name dic...')
        df = pd.read_excel(self.data_folder + 'genes_list_hg19.xlsx')
        dic = {}
        name_set = set()
        for _, row in df.iterrows():
            id_ = row['gene_name']
            name = row['gene_symbol'][1:-1]
            if name in name_set:
                print(name)
                continue
            name_set.add(name)
            dic[id_] = name
            dic[name] = id_
        pkl.dump(dic, open(self.data_folder + 'dump/gene_name.pkl', 'wb'))
        print('done')

    def get_gene_name(self, gene):
        try:
            return self.dic[gene]
        except:
            print(f'{gene} not found')
            return gene


if __name__ == '__main__':
    dic = GeneNameDic('../../data/')
    dic.generate_dic()
