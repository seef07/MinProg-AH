from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms.reinforcement.deepqleaning import run_protein_folding


def main():
	create_experiment_folders()
	sequence = 'HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP'

	best_state, algorithm = run_protein_folding(sequence, 50000)

	algorithm.plot(f'./output/_len.png')
	best_state.plot(f'./output/best_len.png')


if __name__ == '__main__':
	main()
