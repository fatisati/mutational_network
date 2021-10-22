import os

import matplotlib.pyplot as plt

import community_analysis
import config
import data_loader
import plotly_template


class Figures:
    # def __init__(self):
    # self.coms = data_loader.load_coms(0.15)
    # self.patient_gene = data_loader.load_patient_gene()
    # self.patients = list(self.patient_gene.keys())

    def count_covering_coms(self, patient):
        cnt = 0
        for com in self.coms:
            if community_analysis.com_cover_patient(com, self.patient_gene[patient]):
                cnt += 1
        return cnt

    def count_covering_patients(self, com):
        cnt = 0
        for patient in self.patients:
            if community_analysis.com_cover_patient(com, self.patient_gene[patient]):
                cnt += 1
        return cnt

    def sample_coverages(self):
        covered_patient_cnts = []
        for com in self.coms:
            covered_patient_cnts.append(self.count_covering_patients(com))

        com_sizes = [len(com) for com in self.coms]
        covered_com_cnt = [self.count_covering_coms(patient) for patient in self.patients]

        plt.title('number of patients covered by each community')
        plt.xlabel('community size')
        plt.ylabel('number of covered patients')
        plt.plot(com_sizes, covered_patient_cnts)
        plt.show()

        plt.plot(covered_com_cnt)
        plt.show()

    def assortativity(self):
        df = data_loader.load_assortativity_table()
        return plotly_template.show_table(df)

    def degree_dist(self):
        thrs = [0.1, 0.15, 0.2]
        figs = []
        for thr in thrs:
            net = data_loader.load_network(thr)
            degrees = [net.degree(node) for node in net.nodes()]
            figs.append(plotly_template.histogram(degrees))
        return figs

    def go_process(self):
        pass

    def go_function(self):
        pass

    def go_component(self):
        pass

    def cancer_pathways_clustermap(self):
        pass

    def gene_type_piechart(self):
        pass

    def ppi_vs_com_links(self):
        pass

    def coms_ppi_enrichment(self):
        pass

    def cis_trans(self):
        pass


if __name__ == '__main__':
    fig = Figures()
    # print(os.listdir(config.old_table_path))
    table = fig.assortativity()
    plotly_template.show_go(table)