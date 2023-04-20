import pandas as pd
import requests  ## python -m pip install requests
import json
import pickle as pkl
import shared_names
from biological_analysis.string_names import StringNames


class StringApi:
    def __init__(self, res_path):
        string_api_url = "https://version-11-5.string-db.org/api"
        output_format = "json"
        method = "enrichment"
        self.request_url = "/".join([string_api_url, output_format, method])
        self.gene_name_dict = pkl.load(open(res_path + shared_names.gene_name_dic, 'rb'))
        self.res_path = res_path

    def get_enrichment(self, genes, max_fdr=0.05):
        # print(f'getting string enrichment for gene-set with size:{len(genes)}...')

        params = {
            "identifiers": "%0d".join(genes),  # your protein
            "species": 9606,  # species NCBI identifier
            "caller_identity": "www.awesome_app.org"  # your app name
        }
        response = requests.post(self.request_url, data=params)
        data = json.loads(response.text)

        try:
            df = pd.DataFrame(data)
            df['fdr'] = [float(fdr) for fdr in df['fdr']]
            df = df[df['fdr'] <= max_fdr]
            go = df[df['category'].isin(['Component', 'Process', 'Function'])]
            kegg = df[df['category'] == 'KEGG']

            # print('done')
            return list(kegg['term']), list(go['term'])
        except:
            # print(data)
            return [], []

    def coms_gene_enrichment(self, coms):
        go_list, kegg_list = [], []
        for com in coms:
            gene_names = [self.gene_name_dict[gene_id] for gene_id in com]
            go, kegg = self.get_enrichment(gene_names)
            go_list.append(go)
            kegg_list.append(kegg)
        go_cnt_list = [len(go) for go in go_list]
        kegg_cnt_list = [len(kegg) for kegg in kegg_list]
        df = pd.DataFrame({'com': coms, 'go': go_list, 'go-cnt': go_cnt_list,
                           'kegg': kegg_list, 'kegg-cnt': kegg_cnt_list})
        df.to_csv(self.res_path + 'gene-set-enrichment.csv')
        return go_list, kegg_list


if __name__ == '__main__':
    res_path = '../results/'
    coms = pkl.load(open(res_path + 'h-ensemble1000-coms0.15.pkl', 'rb'))
    string_api = StringApi(res_path)
    string_api.coms_gene_enrichment(coms)