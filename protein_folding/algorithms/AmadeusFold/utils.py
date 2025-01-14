import random
import torch
import numpy as np


_molecule_indices = {'H': 0, 'C': 1, 'P': 2}

def fold_validity_quantified_loss(coordinates) -> float:
	"""
	This function returns a quantified measure of how "invalid" a given
	ordering is, multiplied by adjustable weights. This allows invalid
	orderings to still form a gradient, allowing models to move in the right
	direction. 

	Pre:
		- coordinates is a torch tensor containing node coordinates
	Post:
		- a float containing the accumulated loss is returned
	"""
	multiplier_duplicate_coordinates = 10
	multiplier_invalid_steps = 10

	# penalty for duplicate coordinates, ie overlapping molecules
	loss = (multiplier_duplicate_coordinates *
		(len(coordinates[0]) - len(coordinates.unique(dim=1)[0])))
	
	# penalty for sequential molecules that are not spaced 1 unit apart
	for idx in range(1, len(coordinates[0])):
		coor_1 = coordinates[0][idx - 1]
		coor_2 = coordinates[0][idx]
		diff = coor_1 - coor_2
		if not (torch.count_nonzero(diff) == 1 and
				torch.abs(torch.sum(diff)) == 1):
			loss += multiplier_invalid_steps
			print(coor_1, coor_2, diff)

	return loss


def _get_bond_score(letter1: str, letter2: str):
	if letter1 == 'P' or letter2 == 'P':
		return 0
	if letter1 == 'H' and letter2 == 'H':
		return -1
	elif letter1 == 'C' and letter2 == 'C':
		return -5
	else:
		# C-H or H-C, only options left so avoid more comparisons for speed
		return -1


def compute_bond_score_coordinates(seq, coordinates):
	"""
	This variation of the bond score computing function takes the coordinates
	of nodes directly as input, instead of as a list of orders.

	Pre:
		- coordinates is a torch tensor containing node coordinates
	Post:
		- a bond score is returned as float
	"""
	# Keep track of protein score
	score = 0

	# Create dictionary mapping node positions to their letter and add first
	# node
	nodes = {(0,0,0): seq[0]}

	seq_idx = 1
	order_idx = 0
	for coordinate in coordinates[0][1:]:
		x = int(coordinate[0])
		y = int(coordinate[1])
		z = int(coordinate[2])

		# Do not save node or check score if it's a 'P', since it will not
		# affect the final score
		if seq[seq_idx] != 'P':
			nodes[(x, y, z)] = seq[seq_idx]

			# Remove score to compensate for over-counting in advance
			score -= _get_bond_score(seq[seq_idx], seq[seq_idx - 1])

		seq_idx += 1
		order_idx += 1

	# Go through all values in dict and add score if there is a node right, up
	# or forward
	for pos, letter in nodes.items():
		# Check right
		try_pos = (pos[0] + 1, pos[1], pos[2])
		neighbour_letter = nodes.get(try_pos)
		if neighbour_letter:
			score += _get_bond_score(letter, neighbour_letter)

		# Check up
		try_pos = (pos[0], pos[1] + 1, pos[2])
		neighbour_letter = nodes.get(try_pos)
		if neighbour_letter:
			score += _get_bond_score(letter, neighbour_letter)

		# Check forward
		try_pos = (pos[0], pos[1], pos[2] + 1)
		neighbour_letter = nodes.get(try_pos)
		if neighbour_letter:
			score += _get_bond_score(letter, neighbour_letter)

	return score


def index_sequence(sequence):
	return [[_molecule_indices[letter] for letter in sequence]]


def sequence_generator(length: int, fractions: list = [.5, .25, .25]) -> list:
	"""
	generate a sequence of len length of valid molecules (letters). Fractions
	denotes the fraction of the total letters that will be (not exactly)
	composed of that letter. Ex: for fractions = [.5, .25, .25], half of the
	letters in the sequence will be H, a quarter will be C and a quarter will
	be P.
	
	pre:
		- length is an int > 0
		- fractions is a list of len = 3 containing floats, the sum of which
		is 1.0
	post:
		- returns a list of length = length made up of H, C, P
	"""
	sequence = ['H'] * int(length * .5)
	sequence += ['C'] * int(length * .25)
	sequence += ['P'] * (length - int(length * .5) - int(length * .25))
	random.shuffle(sequence)

	return sequence