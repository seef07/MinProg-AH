from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import IterativeRandom


def main():
    create_experiment_folders()

    sequence = 'H' * 30
    sequence = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    protein = Protein(sequence)

    algorithm = IterativeRandom(protein, dimensions=2, debug=True)
    score = algorithm.run()

    protein.plot(f'./output/{algorithm.get_name()}_protein_len{len(sequence)}.png')

    print(f'Final score: {score}')


if __name__ == '__main__':
    main()
