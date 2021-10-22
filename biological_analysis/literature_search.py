import data_loader


def find_row(df, key, val):
    return df[df[key] == val]


if __name__ == '__main__':
    gene = 'MAGI2'
    cancer_df, breast_cancer_df = data_loader.load_pmids()
    print(find_row(cancer_df, 'gene', gene), find_row(breast_cancer_df, 'gene', gene))
