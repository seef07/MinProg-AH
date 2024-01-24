# Import necessary modules
from protein_folding.protein import Protein
from protein_folding.node import Node

# Example Heuristic 1: Distance Heuristic
def distance_heuristic(state):
    # Evaluate state based on the total or average distance between nodes
    total_distance = 0
    for i in range(len(state.nodes) - 1):
        total_distance += state.nodes[i].distance_to(state.nodes[i + 1])
    return total_distance

########################letter p ####################
def currentletter(protein_sequence, fold_index):
    """
    Apply the folding heuristic for a given fold index in the protein sequence.
    """
    # Ensure the fold_index is within valid range
    if fold_index < 2 or fold_index > len(protein_sequence) - 3:
        return 0

    # Extract the pattern around the fold point
    p = protein_sequence[fold_index]

    if p == 'P':
        return 1

##############Energy heuristic#####################################


##############################Compactness Heuristic############
def compactness_heuristic(current_state, next_state):
    def calculate_compactness(state):
        if not state:
            return 0
        avg_x = sum(node[0] for node in state) / len(state)
        avg_y = sum(node[1] for node in state) / len(state)
        variance = sum((node[0] - avg_x)**2 + (node[1] - avg_y)**2 for node in state) / len(state)
        return 1 / (1 + variance)

    return calculate_compactness(next_state) - calculate_compactness(current_state)

###################pattern heuristic#####################
def score_pattern(pattern):
    """
    Score the given pattern based on predefined rules.
    """
    if pattern in ['HH', 'HHHH']:
        return 2
    elif pattern in ['PP', 'PPPP']:
        return 1
    elif pattern in ['HP', 'PH', 'HPHP', 'PHHP']:
        return 0
    elif pattern in ['HHPP', 'PPHH']:
        return -1
    elif pattern in ['HHHH', 'PPPP']:
        return -2
    elif pattern in ['HPHH', 'HHPH']:
        return 1.5
    elif pattern == 'HPHP':
        return -1.5
    else:
        return 0

def folding_heuristic(protein_sequence, fold_index):
    """
    Apply the folding heuristic for a given fold index in the protein sequence.
    """
    # Ensure the fold_index is within valid range
    if fold_index < 2 or fold_index > len(protein_sequence) - 3:
        return 0

    # Extract the pattern around the fold point
    pattern = protein_sequence[fold_index-2 : fold_index+3]

    # Score the pattern
    score = score_pattern(pattern)

    return score
###########################################################