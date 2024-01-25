
###########alt heuristic############
def altheuristic(protein_sequence, fold_index):
    hydrophobic = {'H'}  # Set of hydrophobic residue types
    polar = {'P'}  # Set of polar residue types

    # Define the range to consider around the folding point
    window_size = 4  # Adjust as needed
    start = max(0, fold_index - window_size)
    end = min(len(protein_sequence), fold_index + window_size + 1)

    segment = protein_sequence[start:end]

    hydrophobic_count = 0
    polar_count = 0
    for residue in segment:
        if residue in hydrophobic:
            hydrophobic_count += 1
        elif residue in polar:
            polar_count += 1

    if hydrophobic_count + polar_count == 0:
        return 0  # Avoid division by zero

    score = hydrophobic_count / (hydrophobic_count + polar_count)
    return score

########################letter p ####################
def currentletter(protein_sequence, fold_index):
    """
    Apply the folding heuristic for a given fold index in the protein sequence.
    Returns 1 if the amino acid at the fold index is 'P', otherwise returns 0.
    """
    # Ensure the fold_index is within the valid range
    if fold_index < 2 or fold_index > len(protein_sequence) - 3:
        return 0

    # Extract the amino acid at the fold point
    p = protein_sequence[fold_index]

    # Return 1 if the amino acid is 'P', else return 0
    return 1 if p == 'P' else 0


##############Energy heuristic#####################################


##############################Compactness Heuristic############
def compactness_heuristic(current_state, next_state):
    def calculate_compactness(state):
        if not state:
            return 0
        # Calculate the average coordinates
        avg_x = sum(node[0] for node in state) / len(state)
        avg_y = sum(node[1] for node in state) / len(state)

        # Calculate the variance
        variance = sum((node[0] - avg_x)**2 + (node[1] - avg_y)**2 for node in state) / len(state)

        # Return the inverse of variance as a measure of compactness
        return 1 / (1 + variance)

    # Calculate compactness for current and next state and return the difference
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