from utils import *

if __name__ == '__main__':
    proteins = ['insulin', 'keratin', 'pikachurin']
    prepare_data(proteins)
    data = load_data(get_proteins_filename(proteins, output=True))
    cluster = cluster_data(data, nclusters=len(proteins))
    print(cluster)
    draw_tree(data)
