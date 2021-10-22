import data_loader


class EnrichmentAnalysis:
    def __init__(self):
        # self.go = data_loader.load_coms_ontology()
        self.cancer_pathway_genes = data_loader.load_cancer_pathways_genes()
