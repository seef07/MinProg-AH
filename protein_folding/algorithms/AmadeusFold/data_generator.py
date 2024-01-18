import sys
sys.path.append('../../..')

import pandas as pd
import random

from protein_folding.algorithms.bruteforce import BruteForce
from protein_folding.protein import Protein
from sequence_generator import sequence_generator

def data_generator(
	n_sequences: int = 3,
	seq_len_limits: list = [5, 10],
	max_bruteforce_iterations: int = 10000,
	csv_loc: str = './data/test_data.csv'
	) -> None:

	dataset = list()

	# loop over number of desired sequences
	for _ in range(n_sequences):
		# generate a random sequence of random length within bounds
		n_molecules = random.randint(seq_len_limits[0], seq_len_limits[1])
		sequence = sequence_generator(n_molecules)
		protein = Protein(sequence)

		# generate brute force solutions
		algorithm = BruteForce(protein, dimensions=2, max_iterations=max_bruteforce_iterations)
		results = algorithm.run()

		# obtain data from solutions: sequence, score, order, node coordinates
		for	order in results:
			score = results[order]
			protein.set_order(order)
			coordinates = protein.get_node_coordinates()

			datapoint = [sequence, score, order, coordinates]
			dataset.append(datapoint)
	
	# convert to dataframe and save
	df = pd.DataFrame(dataset, columns=['sequence', 'score', 'order', 'coordinates'])
	df.to_csv(csv_loc)

if __name__ == "__main__":
	data_generator()