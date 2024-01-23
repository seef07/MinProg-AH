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


