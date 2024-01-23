import numpy as np
import matplotlib.pyplot as plt
from numba import jit
import random
import matplotlib
matplotlib.rcParams.update({'font.size': 15})

def define_sequence(seq_str):
    amino_acid_map = {'H': 1, 'P': 2, 'C': 3}
    return [amino_acid_map[acid] for acid in seq_str]

