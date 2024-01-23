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



