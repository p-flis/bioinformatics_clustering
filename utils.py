from Bio.Align import MultipleSeqAlignment
from typing import Dict, List
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import Bio.Phylo
import Bio.Cluster
import Bio.SeqIO
import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess


def get_proteins_filename(proteins: List[str], output: bool = False) -> str:
    """ Returns proteins filename. """
    result = '_'.join(sorted(proteins))
    if output:
        result += '_output'
    return result + '.fasta'


def prepare_data(proteins: List[str]) -> None:
    """ Takes FASTA proteins and merges them to one file. """

    proteins_data = []

    for protein in proteins:
        protein_path = os.path.join("data", protein)
        sequences = os.listdir(protein_path)

        for sequence in sequences:
            with open(os.path.join(protein_path, sequence)) as f:
                for line in Bio.SeqIO.parse(f, "fasta"):
                    line.description = protein + ' ' + line.description
                    line.id = protein + '_' + line.id
                    proteins_data.append(line)

    with open(os.path.join(os.path.join("msa", get_proteins_filename(proteins))), 'w+') as f:
        for protein_data in proteins_data:
            Bio.SeqIO.write(protein_data, f, "fasta")


def run_muscle(protein: str) -> None:
    """ Runs muscle program on selected proteins and saves result to file. """

    subprocess.run(f"muscle.exe -align"
                   f"{os.path.join('msa', protein)}"
                   f"-output "
                   f"{os.path.join('outputs', protein + '_output.fasta')}")


def load_data(protein: str) -> MultipleSeqAlignment:
    """ Loads multiple segment alignment from file (result of run_muscle function)."""

    msa = MultipleSeqAlignment([])
    with open(os.path.join('msa', protein)) as f:
        for line in Bio.SeqIO.parse(f, "fasta"):
            msa.append(line)
    msa.sort(key=lambda x: x.description.split(' ')[0])
    return msa


def cluster_data(msa: MultipleSeqAlignment, nclusters: int = 8) -> Dict[str, int]:
    """ Clusters multiple segment alignment proteins using Bio.Cluster. """

    calculator = DistanceCalculator('identity')
    diffs = calculator.get_distance(msa)
    distance_matrix = np.array(diffs)
    clusters = Bio.Cluster.kcluster(distance_matrix, nclusters=nclusters, npass=10)[0]
    return clusters


def draw_tree(msa: MultipleSeqAlignment) -> None:
    """ Draws tree using prepared multiple segment alignment. """

    constructor = DistanceTreeConstructor()
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(msa)
    tree = constructor.upgma(dm)
    Bio.Phylo.draw(tree)
    plt.show()

