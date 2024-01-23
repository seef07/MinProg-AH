import random
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


# Protein sequence
sequence = "HHPHHHPHPHHHPH"  # Change this to your protein sequence

# Energy matrix
energy_matrix = {
    'HH': -1, 'CC': -5, 'CH': 0, 'HC': 0, 
    'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
}

# Initialize positions in a straight line
positions = [(i, 0) for i in range(len(sequence))]

def calculate_energy(positions, sequence):
    energy = 0
    for i in range(len(sequence)):
        for j in range(i + 2, len(sequence)):  # Start from i+2 to skip consecutive amino acids
            if abs(positions[i][0] - positions[j][0]) + abs(positions[i][1] - positions[j][1]) == 1:
                pair = ''.join(sorted([sequence[i], sequence[j]]))
                energy += energy_matrix[pair]
    return energy


def is_valid_configuration(positions):
    return len(positions) == len(set(positions))

