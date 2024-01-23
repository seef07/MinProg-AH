import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import random
import matplotlib
matplotlib.rcParams.update({'font.size': 15})

def define_sequence(seq_str):
    amino_acid_map = {'H': 1, 'P': 2, 'C': 3}
    return [amino_acid_map[acid] for acid in seq_str]

def initialize(N_aminoacids, seq):
    A = np.array(define_sequence(seq))
    J = np.zeros((4, 4))
    J[1, 1] = -1  # HH interactie
    J[3, 3] = -5  # CC interactie
    return A, J


def step_2D_jit():
    d = np.random.randint(1, 5)
    if d == 1:
        x_s, y_s = 1, 0
    elif d == 2:
        x_s, y_s = 0, 1
    elif d == 3:
        x_s, y_s = -1, 0
    elif d == 4:
        x_s, y_s = 0, -1
    return x_s, y_s

def neighbor(x, y, x_test, y_test):
    return (x - x_test) ** 2 + (y - y_test) ** 2 < 1.1




def Calculate_delta_Energy_jit(N_aminoacids, R, J, A, index, x_s, y_s, temp):
    k = index
    x_s, y_s = x_s, y_s
    R_new = np.copy(R)
    R_new[k:, 0] += x_s  # Verplaats alle x-coördinaten vanaf index k
    R_new[k:, 1] += y_s  # Verplaats alle y-coördinaten vanaf index k

    # Controleer of alle nieuwe posities uniek zijn
    for i in range(len(R_new)):
        for j in range(i + 1, len(R_new)):
            if np.all(R_new[i] == R_new[j]):
                return R 
    
    energy_before = 0
    energy_new = 0

    for i in range(N_aminoacids):
        for j in range(i):
            if i - j > 1:
                distance_before = np.linalg.norm(R[i] - R[j])
                distance_new = np.linalg.norm(R_new[i] - R_new[j])
                if distance_before < 1.1:
                    energy_before -= J[A[i], A[j]]
                if distance_new < 1.1:
                    energy_new -= J[A[i], A[j]]

    delta_energy = energy_new - energy_before
    if delta_energy <= 0 or random.random() < np.exp(-delta_energy / temp):
        R[k:, 0] += x_s
        R[k:, 1] += y_s

    return R


def is_valid_move(R, x_new, y_new, k):
    if k > 0 and np.abs(R[k - 1, 0] - x_new) + np.abs(R[k - 1, 1] - y_new) != 1:
        return False
    for i in range(k):
        if R[i, 0] == x_new and R[i, 1] == y_new:
            return False
    return True

def calculate_new_positions(R, k, x_s, y_s):
    N_aminoacids = len(R)
    new_positions = np.copy(R)

    if k > 0:
        new_positions[k, 0] = new_positions[k-1, 0] + x_s
        new_positions[k, 1] = new_positions[k-1, 1] + y_s

        for i in range(k + 1, N_aminoacids):
            new_positions[i, 0] = new_positions[i-1, 0] + x_s
            new_positions[i, 1] = new_positions[i-1, 1] + y_s

    return new_positions

def calculate_energy(R, J, A):
    energy = 0
    N_aminoacids = len(R)
    for i in range(N_aminoacids):
        for j in range(i):
            if i - j > 1:
                distance = np.linalg.norm(R[i] - R[j])
                if distance < 1.1:
                    energy -= J[A[i], A[j]]
    return energy


def ProteinFolding_jit(N_aminoacids, J, A, Temp=10, Nsteps=1000000):
    R = np.zeros((N_aminoacids, 2))
    for i in range(N_aminoacids):
        R[i, 0], R[i, 1] = i, 0

    FreeEnergy = []
    length = []
    step = []

    for i in range(Nsteps):
        k = np.random.randint(1, N_aminoacids)
        x_s, y_s = np.random.randint(-1, 2), np.random.randint(-1, 2)

        new_positions = calculate_new_positions(R, k, x_s, y_s)
        # Pass x_new and y_new as separate arguments instead of new_positions[k]
        if not is_valid_move(R, new_positions[k, 0], new_positions[k, 1], k):
            print("invalid")
            continue
        F_new = calculate_energy(new_positions, J, A)
        F_old = calculate_energy(R, J, A)
        delta_energy = F_new - F_old

        if delta_energy <= 0 or random.random() < np.exp(-delta_energy / Temp):
            R = new_positions
            visualize_positions(R, A)

        r2 = np.sum((R[N_aminoacids - 1] - R[0]) ** 2)
        step.append(i + 1)
        FreeEnergy.append(F_new)
        length.append(np.sqrt(r2))

    return step, FreeEnergy, length, R


def initialize_Uniform_J(N_aminoacids):
    # Deze functie zal nu uniforme interactie-energieën instellen
    A = np.zeros(N_aminoacids, int)
    J = np.zeros((4, 4)) + 5
    # Hier kun je de waarden van A instellen zoals nodig is voor je toepassing
    return A, J

 ##########################PloT#################################################################           
def plot_folding_structure(R, A):
    """Plot de structuur van het gevouwen eiwit."""
    for i in range(len(R) - 1):
        plt.plot([R[i][0], R[i+1][0]], [R[i][1], R[i+1][1]], 'ro-')
        plt.text(R[i][0], R[i][1], str(A[i]))

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Eiwit Vouwstructuur')
    plt.grid(True)
    plt.show()

def visualize_positions(R, A):
    # Create a map of amino acids for visualization
    amino_acid_symbols = {1: 'H', 2: 'P', 3: 'C'}

    # Determine the size of the grid
    x_min, y_min = np.min(R, axis=0)
    x_max, y_max = np.max(R, axis=0)

    grid_width = int(x_max - x_min + 1)
    grid_height = int(y_max - y_min + 1)

    # Create an empty grid
    grid = [[' ' for _ in range(grid_width)] for _ in range(grid_height)]

    # Place the amino acids on the grid
    for pos, acid in zip(R, A):
        x, y = pos
        grid[int(y - y_min)][int(x - x_min)] = amino_acid_symbols[acid]

    # Print the grid
    for row in grid[::-1]:  # reverse for correct y-axis orientation
        print(' '.join(row))

##############################################################################################
