from protein_folding.protein import Protein
from A_new_datastru.montecarlo import monte_carlo_folding, visualize_protein, plot_weights_and_bond_scores
from A_new_datastru.D3_MonteCarlo import monte_carlo_folding_3, visualize_protein_3

import os
import sys
from protein_folding.protein import Protein
from protein_folding.algorithms import *
from protein_folding.algorithms.heuristics import *


from matplotlib.animation import PillowWriter
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

energy_matrix = {
    'HH': -1, 'CC': -5, 'CH': -1, 'HC': -1, 
    'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
}

def get_sequence() -> str:
    match len(sys.argv):
        case 1:
            sequence = input("Enter sequence: ")
        case 2:
            sequence = sys.argv[2]
        case _:
            print("Usage: python main.py [sequence] or python main.py and input sequence manually.")
            sys.exit(1)

    return sequence

def setup_visualization(positions, sequence, energy_matrix):
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    colors = {'H': 'red', 'P': 'blue', 'C': 'green'}
    color_sequence = [colors[acid] for acid in sequence]

    x, y, z = zip(*positions)
    ax.scatter(x, y, z, c=color_sequence, depthshade=False, s=250)

    for i in range(len(positions) - 1):
        ax.plot([positions[i][0], positions[i + 1][0]], [positions[i][1], positions[i + 1][1]],
                [positions[i][2], positions[i + 1][2]], color='black', linestyle='-', linewidth=4)

    for i in range(len(sequence)):
        for j in range(i + 2, len(sequence)):
            dx = abs(positions[i][0] - positions[j][0])
            dy = abs(positions[i][1] - positions[j][1])
            dz = abs(positions[i][2] - positions[j][2])

            if (dx == 1 and dy == 0 and dz == 0) or (dx == 0 and dy == 1 and dz == 0) or (dx == 0 and dy == 0 and dz == 1):
                pair = ''.join(sorted([sequence[i], sequence[j]]))
                line_color = 'red'  # Default color for interactions
                linestyle = '--'
                if pair == 'CC':  # Specific case for CC pairs
                    line_color = 'green'
                    linestyle = '--'  # Dashed line for CC pairs
                if pair == 'CH' or pair == 'HC':  # Specific case for CC pairs
                    line_color = 'blue'
                    linestyle = '--'  # Dashed line for CC pairs
                if energy_matrix.get(pair, 0) != 0:
                    ax.plot([positions[i][0], positions[j][0]], [positions[i][1], positions[j][1]],
                            [positions[i][2], positions[j][2]], color=line_color, linestyle=linestyle, linewidth=2)
    ax.set_axis_off()

    return fig, ax

def update(num, ax, angle_step):
    ax.view_init(elev=20, azim=num * angle_step)

def create_and_save_animation(positions, sequence, energy_matrix, filename='protein_folding.gif'):
    fig, ax = setup_visualization(positions, sequence, energy_matrix)
    angle_step = 2 
    frames = int(360 / angle_step) 

    animation = FuncAnimation(fig, update, frames=frames, fargs=(ax, angle_step), blit=False)
    writer = PillowWriter(fps=20)  
    animation.save(filename, writer=writer)
    plt.close(fig)  


def main():
    print("Enter before we start a sequence.")
    sequence = get_sequence()
    protein = Protein(sequence)
    print("Protein Folding Simulator")
    print("=========================")
    print("Available Algorithms:")
    print("1. Monte Carlo")
    print("2. Monte Carlo 3D")
    print("3. Bruteforce")
    print("4. Greedy")
    print("5. Iterative Greedy")
    print("6. Iterative Random")
    print("7. Simulated Annealing")
    a_i = input("Enter Algorithm InteGer: ")

    if int(a_i) == 1:
        #MonteCarlo
        best_positions, best_energy = monte_carlo_folding(sequence)
        print("Best positions:", best_positions)
        print("Best energy:", best_energy)
        visualize_protein(best_positions, sequence)
        plot_weights_and_bond_scores()
    elif int(a_i) == 2:
        #montecarlo 3d
        best_positions, best_energy = monte_carlo_folding_3(sequence)
        print("Best energy:", best_energy)
        visualize_protein_3(best_positions, sequence, energy_matrix)
    elif int(a_i) == 3:
    #Bruteforce
        max_iterations = int(input("Choose max iteration: "))
        dim = int(input("Choose dimensions [2/3]: "))
        protein = Protein(sequence)
        algorithm = BruteForce(protein, dimensions=dim, max_iterations=max_iterations, verbose=True)
        results = algorithm.run()
        print(results)
        print(f"Lowest score: {min(results.values())}")
    elif int(a_i) == 4:
        ## Greedy
        dim = int(input("Choose dimensions [2/3]: "))
        protein = Protein(sequence)

        algorithm = Greedy(protein, dimensions=dim, debug=True)
        score = algorithm.run()
        print(f"Score: {score}")
        #if dim == 2:
        #    protein.plot(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
        #elif dim == 3:
        #    protein.plot_3d(f'./output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
    elif int(a_i) == 5:
    #iterative greedy
        max_iterations = int(input("Choose max iteration: "))
        dim = int(input("Choose dimensions [2/3]: "))
        algorithm = IterativeGreedy(protein, dimensions=dim, max_iterations=max_iterations, debug=True)
        score = algorithm.run()
        algorithm.plot_bond_scores()
        print(f"Sequence: {sequence}")
        print(f"Score: {score}")
        #if dim == 2:
        #    protein.plot(f'./experiments/output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
        #elif dim == 3:
        #   protein.plot_3d(f'./experiments/output/evaluate_{algorithm.get_name()}_len{len(sequence)}_dim{dim}.png')
        ########################
    elif int(a_i) == 6:
        #iterativerandom
        dim = int(input("Choose dimensions [2/3]: "))
        algorithm = IterativeRandom(protein, dimensions=dim, debug=True)
        score = algorithm.run()
        print(f'Final score: {score}')

    elif int(a_i) == 7: 
    #simulated annealing
        dim = int(input("Choose dimensions [2/3]: "))
        algorithm = SimulatedAnnealing(protein, dimensions=dim, debug=True)
        score = algorithm.run()
        print(f"Score: {score}")

main()
