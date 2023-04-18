import plotly_template
from visualization import Utils
import pickle as pkl
import shared_names

import pandas as pd


def com_cover_patient(com, patient_genes):
    if len(set(com).intersection(set(patient_genes))) == len(com):
        return True
    return False


class CommunityAnalysis:
    def __init__(self, coms, result_path):
        self.coms = coms
        self.result_path = result_path
        patient_gene_path = f'{result_path}/{shared_names.patient_gene}'
        self.patient_gene = pkl.load(open(patient_gene_path, 'rb'))
        self.patients = list(self.patient_gene.keys())

        self.com_patient_coverage = None
        self.covered_coms_cnt = None

        self.vis = Utils(result_path)

    def covering_coms(self, patient):
        coms = []
        for com in self.coms:
            if com_cover_patient(com, self.patient_gene[patient]):
                coms.append(com)
        return coms

    def count_covering_patients(self, com):
        cnt = 0
        for patient in self.patients:
            if com_cover_patient(com, self.patient_gene[patient]):
                cnt += 1
        return cnt

    def coms_report(self):
        print(f'number of communities: {len(self.coms)}')

        results = {'coms': self.coms}

        print('count covering patients...')
        com_patient_cnt = [self.count_covering_patients(com) for com in self.coms]
        results['covering-patients-cnt'] = com_patient_cnt
        print('done')

        com_patient_coverage = [cnt / len(self.patients) for cnt in com_patient_cnt]
        results['patients-coverage'] = com_patient_coverage

        pd.DataFrame(results).to_csv(self.result_path + 'coms-patient-coverage.csv')
        self.com_patient_coverage = com_patient_coverage

    def patient_report(self):
        print('count patient-com')
        covered_coms = [self.covering_coms(patient) for patient in self.patients]
        self.covered_coms_cnt = [len(coms) for coms in covered_coms]
        results = {'patient': self.patients, 'covered-com-cnt': self.covered_coms_cnt, 'covered-coms': covered_coms}
        not_covered_patients = self.covered_coms_cnt.count(0)
        print(f'number of patients not covered by any community: {not_covered_patients}, ratio: {not_covered_patients / len(self.patients)}')
        pd.DataFrame(results).to_csv(self.result_path + 'patient-com-covered.csv')

    def report(self):
        self.coms_report()
        self.patient_report()

    def plot(self):
        if not self.com_patient_coverage:
            self.coms_report()
        if not self.covered_coms_cnt:
            self.patient_report()
        coms_size = [len(com) for com in self.coms]
        cp_coverage = plotly_template.scatter(coms_size, self.com_patient_coverage)
        pc_coverage_hist = plotly_template.histogram(self.covered_coms_cnt)

        patient_ids = list(range(len(self.patients)))
        pc_coverage = plotly_template.scatter(patient_ids, self.covered_coms_cnt)

        plot_types = ['Scatter', 'Histogram', 'Scatter']
        xlabels = ['Community size', 'Communities count', 'Patient id']
        ylabels = ['Patient coverage',
                   'Covered patients count', 'Number of communities']
        titles = ['a', 'b', 'c', 'd']
        rows, cols = 3, 1
        fig = plotly_template.make_multiple_plots([cp_coverage, pc_coverage_hist, pc_coverage], rows, cols, plot_types,
                                                  titles, xlabels, ylabels)
        plotly_template.save(fig, self.result_path + 'patient-coverage.svg')
        plotly_template.save(fig, self.result_path + 'patient-coverage.png')


if __name__ == '__main__':
    coms = pkl.load(open('../results/h-ensemble10-coms0.15.pkl', 'rb'))
    coms_analyser = CommunityAnalysis(coms, '../results/')
    coms_analyser.patient_report()
    # coms_analyser.report()
    # coms_analyser.plot()
