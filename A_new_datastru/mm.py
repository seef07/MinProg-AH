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

