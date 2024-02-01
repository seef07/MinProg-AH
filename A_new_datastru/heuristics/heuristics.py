
def altheuristic(protein_sequence, fold_index):
    """
    Evaluates the hydrophobic to polar ratio within a specified window around a fold index in a protein sequence.
    
    Parameters:
    protein_sequence (str): The sequence of amino acids in the protein.
    fold_index (int): The index around which the evaluation is centered.
    
    Returns:
    float: The ratio of hydrophobic to total (hydrophobic + polar) residues in the window.
    """
    hydrophobic = {'H'}  # Set of hydrophobic residue types
    polar = {'P'}  # Set of polar residue types

    window_size = 4
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

def currentletter(protein_sequence, fold_index):
    """
    Determines if the amino acid at the fold index is 'P' or 'H', returning a specified score.
    
    Parameters:
    protein_sequence (str): The sequence of amino acids in the protein.
    fold_index (int): The index of the folding point in the sequence.
    
    Returns:
    int: A score based on the type of amino acid at the fold index.
    """
    if fold_index < 2 or fold_index > len(protein_sequence) - 3:
        return 0

    if protein_sequence[fold_index] == 'P':
        return 4
    elif protein_sequence[fold_index] == 'H':
        return 1
    else:
        return 0

def compactness_heuristic(current_state, next_state):
    """
    Compares the compactness of the current and next states of a protein structure.
    
    Parameters:
    current_state (list of tuples): Coordinates and types of amino acids in the current state.
    next_state (list of tuples): Coordinates and types of amino acids in the next state.
    
    Returns:
    float: The difference in compactness between the next and current states.
    """
    def calculate_compactness(state):
        if not state:
            return 0
        avg_x = sum(node[0] for node in state) / len(state)
        avg_y = sum(node[1] for node in state) / len(state)
        variance = sum((node[0] - avg_x)**2 + (node[1] - avg_y)**2 for node in state) / len(state)
        return 1 / (1 + variance)

    return calculate_compactness(next_state) - calculate_compactness(current_state)


def score_pattern(pattern):
    """
    Score the given protein pattern based on predefined rules. This function is used
    to evaluate a segment of a protein sequence for its structural propensity based on
    the composition and arrangement of its amino acids.

    Parameters:
    - pattern (str): A substring of the protein sequence, typically around a focal point or fold index.

    Returns:
    - int or float: A score representing the evaluated structural propensity of the pattern.
      High positive scores indicate favorable configurations, while negative scores indicate
      less favorable ones.

    The scoring rules are based on hypothetical preferences for certain amino acid arrangements.
    For example, 'HH' or 'HHHH' represent hydrophobic interactions that are favorable in protein
    folding, thus receiving higher scores.
    """
    # Define scoring based on specific patterns
    if pattern in ['HH', 'HHHH']:
        return 4
    elif pattern in ['PP', 'PPPP']:
        return 3
    elif pattern in ['HP', 'PH', 'HPHP', 'PHHP']:
        return 1
    elif pattern in ['HHPP', 'PPHH']:
        return -2
    elif pattern in ['HHHHPPPP', 'PPPPHHHH']:  # Corrected to match function comment examples
        return -3.5
    elif pattern in ['HPHH', 'HHPH']:
        return 2
    elif pattern == 'HPHP':
        return 3
    else:
        return 0

def folding_heuristic(protein_sequence, fold_index):
    """
    Apply a folding heuristic for a specific index within a protein sequence. This heuristic
    evaluates the potential folding propensity based on the amino acid pattern surrounding
    the fold index.

    Parameters:
    - protein_sequence (str): The entire protein sequence being evaluated.
    - fold_index (int): The index within the protein sequence around which the folding propensity
      is to be evaluated.

    Returns:
    - int or float: A score indicating the folding propensity of the specified segment within
      the protein sequence. The score is determined by evaluating the pattern of amino acids
      surrounding the fold index using predefined rules.

    Note:
    The function ensures that the fold_index is within the valid range, avoiding edge cases
    where the fold_index might be too close to the start or end of the protein sequence to
    form a meaningful pattern for evaluation.
    """
    # Ensure the fold_index is within valid range
    if fold_index < 2 or fold_index > len(protein_sequence) - 3:
        return 0

    # Extract the pattern around the fold point and score it
    pattern = protein_sequence[fold_index-2 : fold_index+3]
    score = score_pattern(pattern)

    return score

def distance_heuristic(state, new_state):
    """
    Calculates the Euclidean distance between the first and last coordinates in the state.
    This heuristic can be used to assess the compactness or extension of the protein structure.
    
    :param state: The current state of the protein, represented as a list of tuples (x, y, amino_acid).
    :return: The Euclidean distance between the first and last coordinates.
    """

    if not state or not new_state:
        return 0

    # Extract the first and last coordinates from the state
    first_coordinate = state[0][:2]  # (x, y) of the first element
    last_coordinate = state[-1][:2]  # (x, y) of the last element

    # Calculate the Euclidean distance
    distance = ((last_coordinate[0] - first_coordinate[0])**2 + (last_coordinate[1] - first_coordinate[1])**2)**0.5


        # Extract the first and last coordinates from the state
    f_coordinate = new_state[0][:2]  # (x, y) of the first element
    l_coordinate = new_state[-1][:2]  # (x, y) of the last element

    # Calculate the Euclidean distance
    distance_new = ((l_coordinate[0] - f_coordinate[0])**2 + (l_coordinate[1] - f_coordinate[1])**2)**0.5
    score = distance - distance_new
    return 10*score

def calculate_hydrophobic_compactness(state):
    """
    Calculate the compactness of hydrophobic residues in the protein structure.
    Compactness is measured as the inverse of the average distance of 'H' residues from their centroid.
    """
    # Extract hydrophobic residues
    h_residues = [coord for coord in state if coord[2] == 'H']

    if not h_residues:
        return 0  # Avoid division by zero if there are no hydrophobic residues

    # Calculate centroid of hydrophobic residues
    centroid_x = sum(coord[0] for coord in h_residues) / len(h_residues)
    centroid_y = sum(coord[1] for coord in h_residues) / len(h_residues)

    # Calculate average distance from centroid
    avg_distance = sum(((coord[0] - centroid_x)**2 + (coord[1] - centroid_y)**2)**0.5 for coord in h_residues) / len(h_residues)

    # Compactness is the inverse of average distance
    compactness = 1 / avg_distance if avg_distance != 0 else 0
    return 10*(compactness)

def hydrophobic_compactness_heuristic(current_state, new_state):
    """
    Calculate the difference in compactness of hydrophobic residues between the current and new states.
    """
    current_compactness = calculate_hydrophobic_compactness(current_state)
    new_compactness = calculate_hydrophobic_compactness(new_state)
    return 10*(new_compactness - current_compactness)

def calculate_cytosine_compactness(state):
    """
    Calculate the compactness of cytosine ('C') residues in the protein structure.
    Compactness is measured as the inverse of the average distance of 'C' residues from their centroid.
    """
    # Extract cytosine residues
    c_residues = [coord for coord in state if coord[2] == 'C']

    if not c_residues:
        return 0  # Avoid division by zero if there are no cytosine residues

    # Calculate centroid of cytosine residues
    centroid_x = sum(coord[0] for coord in c_residues) / len(c_residues)
    centroid_y = sum(coord[1] for coord in c_residues) / len(c_residues)

    # Calculate average distance from centroid
    avg_distance = sum(((coord[0] - centroid_x)**2 + (coord[1] - centroid_y)**2)**0.5 for coord in c_residues) / len(c_residues)

    # Compactness is the inverse of average distance
    compactness = 1 / avg_distance if avg_distance != 0 else 0
    return  10* compactness

def cytosine_compactness_heuristic(current_state, new_state):
    """
    Calculate the difference in compactness of cytosine ('C') residues between the current and new states.
    """
    current_compactness = calculate_cytosine_compactness(current_state)
    new_compactness = calculate_cytosine_compactness(new_state)
    return 10 * (new_compactness - current_compactness)

def compactness_heuristic_H(current_state, new_state):
    """
    Calculate a score representing the change in compactness of hydrophobic (H) amino acids
    between the current state and the new state of a protein structure.
    :param current_state: List of tuples representing the current state (x, y, amino_acid).
    :param new_state: List of tuples representing the new state (x, y, amino_acid).
    :return: A single score representing the change in compactness of hydrophobic amino acids.
    """

    def average_distance(state):
        hydrophobic_positions = [node[:2] for node in state if node[2] == 'H']
        if len(hydrophobic_positions) < 2:
            return 0
        
        avg_dist = sum(((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)**0.5 
                       for i, pos1 in enumerate(hydrophobic_positions) 
                       for pos2 in hydrophobic_positions[i+1:]) / len(hydrophobic_positions)
        return avg_dist

    current_avg_dist = average_distance(current_state)
    new_avg_dist = average_distance(new_state)

    # The score is the difference in average distances
    # A positive score indicates increased compactness in the new state
    score = current_avg_dist - new_avg_dist
    return score

def compactness_heuristic_C(current_state, new_state):
    """
    Calculate a score representing the change in compactness of hydrophobic (H) amino acids
    between the current state and the new state of a protein structure.
    :param current_state: List of tuples representing the current state (x, y, amino_acid).
    :param new_state: List of tuples representing the new state (x, y, amino_acid).
    :return: A single score representing the change in compactness of hydrophobic amino acids.
    """

    def average_distance(state):
        hydrophobic_positions = [node[:2] for node in state if node[2] == 'C']
        if len(hydrophobic_positions) < 2:
            return 0
        
        avg_dist = sum(((pos1[0] - pos2[0])**2 + (pos1[1] - pos2[1])**2)**0.5 
                       for i, pos1 in enumerate(hydrophobic_positions) 
                       for pos2 in hydrophobic_positions[i+1:]) / len(hydrophobic_positions)
        return avg_dist

    current_avg_dist = average_distance(current_state)
    new_avg_dist = average_distance(new_state)

    # The score is the difference in average distances
    # A positive score indicates increased compactness in the new state
    score = current_avg_dist - new_avg_dist
    return score