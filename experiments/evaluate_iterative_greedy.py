from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import IterativeGreedy
from random import shuffle
import cProfile
import time

def main():
    create_experiment_folders()
    sequence = 'H' * 21
    sequence += 'C' * 3
    sequence += 'P' * 26
    l = list(sequence)
    shuffle(l)
    sequence = ''.join(l)
    sequence = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"
    score = 0
    dim = 3
    protein = Protein(sequence)
    start_time = time.time()
    algorithm = IterativeGreedy(protein, dimensions=dim, max_iterations=100000, debug=True)
    score = algorithm.run()
    algorithm.plot_bond_scores()
    print(f"Sequence: {sequence}")
    print(f"Score: {score}")
    end_time = time.time()
    if dim == 2:
        protein.plot(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
    elif dim == 3:
        protein.plot_3d(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
    print(f"Elapsed time: {end_time - start_time} seconds")  # Print the elapsed time


if __name__ == '__main__':
    profiler = cProfile.Profile()
    profiler.runcall(main)
    profiler.dump_stats('output/protein.prof')