import matplotlib.pyplot as plt
import itertools

def visualize_protein():
    state = [[0, 0, 'H'], [1, 0, 'C'], [2, 0, 'P'], [3, 0, 'H'], [4, 0, 'P'], [4, -1, 'C'], [4, -2, 'P'], [4, -3, 'H'], [3, -3, 'P'], [2, -3, 'C'], [1, -3, 'H'], [0, -3, 'C'], [-1, -3, 'H'], [-2, -3, 'P'], [-3, -3, 'H'], [-4, -3, 'P'], [-5, -3, 'P'], [-6, -3, 'P'], [-7, -3, 
'H'], [-8, -3, 'P'], [-9, -3, 'P'], [-10, -3, 'P'], [-11, -3, 'H'], [-12, -3, 'P'], [-13, -3, 'P'], [-13, -2, 'P'], [-12, -2, 'P'], [-11, -2, 'H'], [-10, -2, 'P'], [-9, -2, 'C'], [-8, -2, 'P'], [-7, -2, 'H'], [-6, -2, 'P'], [-5, -2, 'P'], [-4, -2, 'P'], [-3, -2, 
'H'], [-2, -2, 'P'], [-1, -2, 'H'], [0, -2, 'H'], [1, -2, 'H'], [2, -2, 'C'], [3, -2, 'C'], [3, -1, 'H'], [2, -1, 'C'], [1, -1, 'H'], [0, -1, 'C'], [-1, -1, 'H'], [-1, 0, 'C'], [-1, 1, 'H'], [0, 1, 'H']]
    # Extract positions and sequence from the state
    positions = [(x, y) for x, y, _ in state]
    sequence = [acid for _, _, acid in state]
    print(positions)
    # Map each amino acid type to a color
    colors = {'H': 'red', 'P': 'blue', 'C': 'green'}
    color_sequence = [colors[acid] for acid in sequence]

    x, y = zip(*positions)
    plt.figure(figsize=(10, 10))
    plt.scatter(x, y, color=color_sequence, s=500)

    # Draw the primary structure connections


    # Function to check if two points are adjacent (but not diagonally)
    def is_adjacent(p1, p2):
        return (p1[0] == p2[0] and abs(p1[1] - p2[1]) == 1) or (p1[1] == p2[1] and abs(p1[0] - p2[0]) == 1)

    # Draw the scoring connections for HH pairs
    for (pos1, acid1), (pos2, acid2) in itertools.combinations(zip(positions, sequence), 2):
        if acid1 == 'H' and acid2 == 'H' and is_adjacent(pos1, pos2):
            plt.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], color='red', linestyle='--', linewidth=4 )
    for (pos1, acid1), (pos2, acid2) in itertools.combinations(zip(positions, sequence), 2):
        if acid1 == 'C' and acid2 == 'C' and is_adjacent(pos1, pos2):
            plt.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], color='green', linestyle='--', linewidth=4 )
    for (pos1, acid1), (pos2, acid2) in itertools.combinations(zip(positions, sequence), 2):
        if acid1 == 'H' and acid2 == 'C' and is_adjacent(pos1, pos2):
            plt.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], color='grey', linestyle='--', linewidth=4 )
    for (pos1, acid1), (pos2, acid2) in itertools.combinations(zip(positions, sequence), 2):
        if acid1 == 'C' and acid2 == 'H' and is_adjacent(pos1, pos2):
            plt.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], color='grey', linestyle='--', linewidth=4 )
    for i in range(len(positions) - 1):
        plt.plot([positions[i][0], positions[i + 1][0]], [positions[i][1], positions[i + 1][1]], color='black', linewidth=5)

    plt.title(f"Protein Folding Visualization\nEnergy: ")
    plt.grid(True)  # Add grid lines
    plt.show()

visualize_protein()