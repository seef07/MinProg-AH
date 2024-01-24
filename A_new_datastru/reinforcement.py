import random
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sequence = "HHPHHHPHPHHHPH"

energy_matrix = {
    ##
}

policy_weights = [0.2,0.2,0.2,0.2,0.2]

episodes = 1000

data = [(i, 0, amino_acid) for i, amino_acid in enumerate(sequence)]

print(data[1])

def dompute_heuristic_score(state, action, policy_weights):
    total = 1
    return total

def reward_function(positions, sequence):
    energy = 0
    for i in range(len(sequence)):
        for j in range(i + 2, len(sequence)):  # Start from i+2 to skip consecutive amino acids
            if abs(positions[i][0] - positions[j][0]) + abs(positions[i][1] - positions[j][1]) == 1:
                pair = ''.join(sorted([sequence[i], sequence[j]]))
                energy += energy_matrix[pair]
    return energy