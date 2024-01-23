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

def rotate_segment(positions, index, clockwise=True):
    if index <= 0 or index >= len(positions) - 1:
        return positions

    pivot = positions[index]
    new_positions = positions.copy()
    for i in range(index + 1, len(positions)):
        dx, dy = positions[i][0] - pivot[0], positions[i][1] - pivot[1]
        if clockwise:
            new_positions[i] = (pivot[0] - dy, pivot[1] + dx)
        else:
            new_positions[i] = (pivot[0] + dy, pivot[1] - dx)

    if is_valid_configuration(new_positions):
        return new_positions
    return positions


def monte_carlo_folding(sequence, iterations=100000, start_temp=1.0, end_temp=0.01):
    current_positions = [(i, 0) for i in range(len(sequence))]
    current_energy = calculate_energy(current_positions, sequence)
    best_positions, best_energy = current_positions, current_energy

    for iteration in range(iterations):
        temp = start_temp - iteration * (start_temp - end_temp) / iterations
        index = random.randint(1, len(sequence) - 2)
        new_positions = rotate_segment(current_positions, index, random.choice([True, False]))
        new_energy = calculate_energy(new_positions, sequence)
        print(new_positions)
        if new_energy < current_energy:
            current_positions, current_energy = new_positions, new_energy
            if new_energy < best_energy:
                best_positions, best_energy = new_positions, new_energy
        elif random.random() < math.exp((current_energy - new_energy) / temp):
            current_positions, current_energy = new_positions, new_energy

    return best_positions, best_energy

def visualize_protein(positions, sequence):
    # Map each amino acid type to a color
    colors = {'H': 'red', 'P': 'blue', 'C': 'green'}
    color_sequence = [colors[acid] for acid in sequence]

    x, y = zip(*positions)
    plt.figure(figsize=(10, 10))
    plt.scatter(x, y, color=color_sequence)

    for i in range(len(positions) - 1):
        plt.plot([positions[i][0], positions[i + 1][0]], [positions[i][1], positions[i + 1][1]], color='black')

    plt.title(f"Protein Folding Visualization\nEnergy: {calculate_energy(positions, sequence)}")
    plt.show()



# Run the simulation
best_positions, best_energy = monte_carlo_folding(sequence)
print("Best positions:", best_positions)
print("Best energy:", best_energy)
