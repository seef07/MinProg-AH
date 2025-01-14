from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import Greedy
from random import shuffle


def main():
    create_experiment_folders()
    sequence = 'H' * 40
    sequence += 'C' * 40
    sequence += 'P' * 20
    l = list(sequence)
    shuffle(l)
    sequence = ''.join(l)
    sequence = "PPPHHPPHHPPPPPHHHHHHHPPHHPPPPHHPPHPP"
    score = 0
    dim = 2
    protein = Protein(sequence)

    algorithm = Greedy(protein, dimensions=dim, debug=True)
    score = algorithm.run()
    print(f"Score: {score}")

    if dim == 2:
        protein.plot(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
    elif dim == 3:
        protein.plot_3d(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')


if __name__ == '__main__':
    main()