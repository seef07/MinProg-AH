import random
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np

sequence = "HHPHPHPHPHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"

energy_matrix = {
    'HH': -1, 'CC': -5, 'CH': -1, 'HC': -1, 
    'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
}

positions = None

def calculate_energy(positions, sequence):
    """
    Calculate the energy of a protein configuration based on interactions between non-adjacent amino acids.
    
    Parameters:
    - positions (list of tuples): Coordinates of amino acids in the protein.
    - sequence (str): Sequence of amino acids in the protein.
    
    Returns:
    - int: The total energy of the configuration.
    """
    energy = 0
    for i in range(len(sequence)):
        for j in range(i + 2, len(sequence)):
            distance = abs(positions[i][0] - positions[j][0]) + \
                       abs(positions[i][1] - positions[j][1]) + \
                       abs(positions[i][2] - positions[j][2])
                       
            if distance == 1:
                pair = ''.join(sorted([sequence[i], sequence[j]]))
                energy += energy_matrix.get(pair, 0)
    return energy

def is_valid_configuration(positions):
    """
    Check if a configuration of positions represents a valid, non-overlapping state.
    
    Parameters:
    - positions (list of tuples): Coordinates of amino acids in the protein.
    
    Returns:
    - bool: True if the configuration is valid, False otherwise.
    """
    return len(positions) == len(set(positions))

def rotate_segment_3d(positions, index, direction, axis):
    """
    Rotate a segment of the protein around a specified axis and direction.
    
    Parameters:
    - positions (list of tuples): Coordinates of amino acids in the protein.
    - index (int): Index of the pivot point for rotation.
    - direction (str): 'clockwise' or 'counterclockwise' rotation.
    - axis (str): Axis of rotation ('x', 'y', or 'z').
    
    Returns:
    - list of tuples: New coordinates after rotation.
    """
    new_positions = positions.copy()
    pivot = positions[index]

    # Rotation logic depending on the axis and direction
    for i in range(index + 1, len(positions)):
        x, y, z = positions[i]
        dx, dy, dz = x - pivot[0], y - pivot[1], z - pivot[2]

        # Apply rotation based on the axis and direction
        if axis == 'x':
            if direction == 'clockwise':
                new_positions[i] = (x, pivot[1] - dz, pivot[2] + dy)
            else:  # counterclockwise
                new_positions[i] = (x, pivot[1] + dz, pivot[2] - dy)
        elif axis == 'y':
            if direction == 'clockwise':
                new_positions[i] = (pivot[0] + dz, y, pivot[2] - dx)
            else:  # counterclockwise
                new_positions[i] = (pivot[0] - dz, y, pivot[2] + dx)
        elif axis == 'z':
            if direction == 'clockwise':
                new_positions[i] = (pivot[0] - dy, pivot[1] + dx, z)
            else:  # counterclockwise
                new_positions[i] = (pivot[0] + dy, pivot[1] - dx, z)

    return new_positions if is_valid_configuration(new_positions) else positions

def monte_carlo_folding_3(sequence, iterations=100000, start_temp=1.0, end_temp=0.01):
    """
    Perform a Monte Carlo simulation to fold a protein into a configuration that minimizes its energy.
    
    Parameters:
    - sequence (str): Sequence of amino acids in the protein.
    - iterations (int): Number of iterations for the simulation.
    - start_temp (float): Starting temperature for simulated annealing.
    - end_temp (float): Ending temperature for simulated annealing.
    
    Returns:
    - tuple: The best positions found and the corresponding energy.
    """
    positions = [(i, 0, 0) for i in range(len(sequence))]
    current_positions = positions
    current_energy = calculate_energy(current_positions, sequence)
    best_positions, best_energy = current_positions, current_energy

    for iteration in range(iterations):
        print(iteration)
        temp = start_temp - (start_temp - end_temp) * (iteration / iterations)
        index = random.randint(1, len(sequence) - 2)
        axis = random.choice(['x', 'y', 'z'])
        direction = random.choice(['clockwise', 'counterclockwise'])
        new_positions = rotate_segment_3d(current_positions, index, direction, axis)
        new_energy = calculate_energy(new_positions, sequence)

        if new_energy < best_energy:
            best_positions, best_energy = new_positions, new_energy

        if new_energy < current_energy or random.random() < math.exp((current_energy - new_energy) / temp):
            current_positions, current_energy = new_positions, new_energy

    return best_positions, best_energy

def visualize_protein_3(positions, sequence, energy_matrix):
    """
    Visualize the 3D structure of the folded protein.
    
    Parameters:
    - positions (list of tuples): Coordinates of amino acids in the protein.
    - sequence (str): Sequence of amino acids in the protein.
    - energy_matrix (dict): Interaction energies between amino acid pairs.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    colors = {'H': 'red', 'P': 'blue', 'C': 'green'}
    color_sequence = [colors[acid] for acid in sequence]

    x, y, z = zip(*positions)
    ax.scatter(x, y, z, c=color_sequence, depthshade=False, s=250)

    for i in range(len(positions)-1):
        ax.plot([positions[i][0], positions[i+1][0]], [positions[i][1], positions[i+1][1]], [positions[i][2], positions[i+1][2]], color='black', linestyle='-', linewidth=5)

    plt.title("Protein Folding Visualization")
    plt.show()
