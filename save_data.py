from utils import *
import pandas as pd
import config


def save_coms(coms, name):
    coms_names = Dictionary().get_coms_names(coms)
    df = []

    for i in range(len(coms_names)):
        for gene in coms_names[i]:
            df.append({'gene': gene['name'], 'community id': i, 'type': gene['type']})
    pd.DataFrame(df).to_excel(config.com_path + name + '.xlsx')
