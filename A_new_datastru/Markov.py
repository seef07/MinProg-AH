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
