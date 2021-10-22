import data_loader

def find_neighbors(gene, network):
    neighbors = list(network[gene].keys())
    return neighbors

if __name__ == '__main__':
    thr = 0.15
    gene = 'MAGI2'

    network = data_loader.load_name_network(thr)
    neighbors = find_neighbors(gene, network)

    print(f'this gene is connected to {len(neighbors)} genes and these genes are {neighbors}')