import networkx as nx
import pickle as pkl


class PpiNetGenerator:
    def __init__(self, data_folder, result_path):
        # https://string-db.org/cgi/download?sessionId=b0rskhtvOySj&species_text=Homo+sapiens
        self.link_file_path = data_folder + '9606.protein.links.v11.5.txt'
        self.protein_info_path = data_folder + '9606.protein.info.v11.5.txt'
        self.gene_dic = self.generate_gene_dic()

        self.result_path = result_path

        self.low_confidence_score = 150
        self.medium_score = 400
        self.high_score = 700

    def get_score_name(self, score):
        if score <= self.low_confidence_score:
            return 'low'
        if score <= self.medium_score:
            return 'medium'
        return 'high'

    def generate_network(self, min_score):
        print(f'generating ppi network with min_score {min_score}...')
        f = open(self.link_file_path)
        line = f.readline()
        g = nx.Graph()
        while line:
            line = f.readline()
            try:
                p1, p2, score = line.split()[:3]
                g1, g2 = self.get_protein_gene(p1), self.get_protein_gene(p2)
                score = float(score)
                if score >= min_score:
                    g.add_edge(g1, g2, weight=score)
            except:
                print(f'cant split line: --{line}--')
        print(f'done. number of nodes: {len(g.nodes)}, number of edges: {len(g.edges)}')
        print('saving dump...')
        pkl.dump(g, open(f'{self.result_path}/ppi-net-{self.get_score_name(min_score)}.pkl', 'wb'))
        print('done')

    def generate_gene_dic(self):
        name_dic = {}
        print('generating protein id-name dic...')
        f = open(self.protein_info_path)
        line = f.readline()
        while line:
            line = f.readline()
            try:
                id_, name = line.split()[:2]
                name_dic[id_] = name
                name_dic[name] = id_
            except:
                print(f'cant split line: --{line}--')
        print(f'done. dic size: {len(name_dic)}')
        return name_dic

    def get_protein_gene(self, protein):
        try:
            return self.gene_dic[protein]
        except:
            print(f'{protein} not found in gene_dic')
            return protein


if __name__ == '__main__':
    netgen = PpiNetGenerator('../../data/', '../../results/')
    netgen.generate_network(400)
