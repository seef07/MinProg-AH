import random
import math
import matplotlib.pyplot as plt

# Define the protein sequence and the energy matrix for interaction between amino acids.
sequence = "HHPHPHPHPHHHHPHPPPHPPPHPPPPHPPPHPPPHPHHHHPHPHPHPHH"  # Example protein sequence

energy_matrix = {
    'HH': -1, 'CC': -5, 'CH': -1, 'HC': -1, 
    'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
}

# Lists to store the history of bond scores during the simulation.
bond_scores_history = []

# Initialize the protein in a straight line configuration.
positions = [(i, 0) for i in range(len(sequence))]

def calculate_energy(positions, sequence):
    """
    Calculate the total energy of a given protein configuration based on its positions and sequence.
    
    Parameters:
    - positions (list of tuples): The coordinates of each amino acid in the protein.
    - sequence (str): The sequence of amino acids in the protein.
    
    Returns:
    - int: The total energy of the configuration.
    """
    energy = 0
    for i in range(len(sequence)):
        for j in range(i + 2, len(sequence)):  # Skip adjacent amino acids to consider non-adjacent interactions only.
            if (positions[i][0] == positions[j][0] and abs(positions[i][1] - positions[j][1]) == 1) or \
               (positions[i][1] == positions[j][1] and abs(positions[i][0] - positions[j][0]) == 1):
                pair = ''.join(sorted([sequence[i], sequence[j]]))
                energy += energy_matrix[pair]
    return energy

def is_valid_configuration(positions):
    """
    Check if a given set of positions represents a valid protein configuration without overlaps.
    
    Parameters:
    - positions (list of tuples): The coordinates of each amino acid in the protein.
    
    Returns:
    - bool: True if the configuration is valid, False otherwise.
    """
    return len(positions) == len(set(positions))

def rotate_segment(positions, index, clockwise=True):
    """
    Rotate a segment of the protein around a pivot point either clockwise or counterclockwise.
    
    Parameters:
    - positions (list of tuples): The current positions of the protein.
    - index (int): The pivot point for rotation.
    - clockwise (bool): Direction of rotation. True for clockwise, False for counterclockwise.
    
    Returns:
    - list of tuples: The new positions of the protein after rotation.
    """
    if index <= 0 or index >= len(positions) - 1:
        return positions  # No rotation possible for first or last position.

    pivot = positions[index]
    new_positions = positions.copy()
    for i in range(index + 1, len(positions)):
        dx, dy = positions[i][0] - pivot[0], positions[i][1] - pivot[1]
        if clockwise:
            new_positions[i] = (pivot[0] - dy, pivot[1] + dx)
        else:
            new_positions[i] = (pivot[0] + dy, pivot[1] - dx)

    return new_positions if is_valid_configuration(new_positions) else positions

def monte_carlo_folding(sequence, iterations=70000, start_temp=0.8, end_temp=0.01):
    """
    Simulate protein folding using the Monte Carlo method.
    
    Parameters:
    - sequence (str): The sequence of amino acids in the protein.
    - iterations (int): The number of iterations to perform.
    - start_temp (float): The starting temperature for the simulation.
    - end_temp (float): The ending temperature for the simulation.
    
    Returns:
    - tuple: The best positions found and the corresponding energy.
    """
    current_positions = [(i, 0) for i in range(len(sequence))]
    current_energy = calculate_energy(current_positions, sequence)
    best_positions, best_energy = current_positions, current_energy

    for iteration in range(iterations):
        temp = start_temp - iteration * (start_temp - end_temp) / iterations
        index = random.randint(1, len(sequence) - 2)
        new_positions = rotate_segment(current_positions, index, random.choice([True, False]))
        new_energy = calculate_energy(new_positions, sequence)
        bond_scores_history.append(new_energy)

        if new_energy < current_energy or random.random() < math.exp((current_energy - new_energy) / temp):
            current_positions, current_energy = new_positions, new_energy
            if new_energy < best_energy:
                best_positions, best_energy = new_positions, new_energy

    return best_positions, best_energy

def visualize_protein(positions, sequence):
    """
    Visualize the protein configuration using matplotlib.
    
    Parameters:
    - positions (list of tuples): The coordinates of each amino acid in the protein.
    - sequence (str): The sequence of amino acids in the protein.
    """
    colors = {'H': 'red', 'P': 'blue', 'C': 'green'}
    color_sequence = [colors[acid] for acid in sequence]

    x, y = zip(*positions)
    plt.figure(figsize=(10, 10))
    plt.scatter(x, y, color=color_sequence)
    for i in range(len(positions) - 1):
        plt.plot([positions[i][0], positions[i + 1][0]], [positions[i][1], positions[i + 1][1]], color='black')
    plt.title("Protein Folding Visualization")
    plt.show()

def plot_weights_and_bond_scores():
    """
    Plot the history of bond scores over the simulation iterations.
    """
    plt.figure(figsize=(12, 6))
    plt.plot(bond_scores_history, label='Bond Score')
    plt.title('Bond Scores Over Iterations')
    plt.xlabel('Iteration')
    plt.ylabel('Bond Score')
    plt.legend()
    plt.show()

# Running the Monte Carlo simulation to fold the protein and visualize the results.
best_positions, best_energy = monte_carlo_folding(sequence)
print("Best positions:", best_positions)
print("Best energy:", best_energy)
visualize_protein(best_positions, sequence)
plot_weights_and_bond_scores()
