from preprocess_data import DataPreprocessor
from network.generator import MutationalNetwork
from network.evaluator import NetworkEvaluator
import community_.detection as community_detection
from community_.evaluator import CommunityAnalysis
from biological_analysis.ppi.ppi_analysis import PpiAnalysis
from biological_analysis.gene_enrichment import StringApi
from biological_analysis.gene_locations import GeneLoc

if __name__ == '__main__':
    data_folder = './data/'
    res_folder = './results/'

    mutation_data = 'simple_somatic_mutation.open.tsv'
    data_preprocessor = DataPreprocessor(data_folder, res_folder, mutation_data)
    data_preprocessor.generate_all()

    threshold = 0.15
    mutnet = MutationalNetwork(res_folder)
    network = mutnet.generate_network(threshold)

    network_evaluator = NetworkEvaluator(network, data_folder)
    network_evaluator.evaluate()

    ensemble_size = 10
    coms = community_detection.hierarchical_ensemble(network, ensemble_size)

    community_analyser = CommunityAnalysis(coms, res_folder)
    community_analyser.report()
    community_analyser.plot()

    ppi_confidence_level = 'medium'
    ppi_analysis = PpiAnalysis(ppi_confidence_level, res_folder)
    ppi_analysis.count_network_ppi_links(network)
    random_cnt = 10
    ppi_analysis.evaluate_coms_ppi_links(network, coms, random_cnt)

    enrichment = StringApi(res_folder)
    enrichment.coms_gene_enrichment(coms)

    gene_location = GeneLoc(res_folder)
    gene_location.analyze_all_coms(network, coms)
