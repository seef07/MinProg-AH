from experiments_helper import create_experiment_folders
from protein_folding.protein import Protein
from protein_folding.algorithms import SimulatedAnnealing
import time

def main():
	create_experiment_folders()
	sequence = 'HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP'
	protein = Protein(sequence)
	start_time = time.time()
	algorithm = SimulatedAnnealing(protein, dimensions=2, debug=True)
	score = algorithm.run()
	print(f"Score: {score}")
	print(f"Elapsed time: {time.time() - start_time} seconds")  # Print the elapsed time

	protein.plot(f'./output/simulated.png')
	

if __name__ == '__main__':
	main()
