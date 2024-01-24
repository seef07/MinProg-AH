import random
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from heuristics import altheuristic, currentletter, compactness_heuristic, folding_heuristic
sequence = "HHPHHHPHPHHHPH"

energy_matrix = {
    ##
}

policy_weights = [0.2,0.2,0.2,0.2,0.2]

episodes = 1000

data = [(i, 0, amino_acid) for i, amino_acid in enumerate(sequence)]


def reward_function(positions, sequence): ##check
    energy = 0
    for i in range(len(sequence)):
        for j in range(i + 2, len(sequence)):  # Start from i+2 to skip consecutive amino acids
            if abs(positions[i][0] - positions[j][0]) + abs(positions[i][1] - positions[j][1]) == 1:
                pair = ''.join(sorted([sequence[i], sequence[j]]))
                energy += energy_matrix[pair]
    return energy

def compute_reward(new_state, current_state): ##check
    positions = []
    sequence = []
    for part in current_state:
        positions.append((part[0], part[1]))
        sequence.append(part[2])

    new_positions = []
    new_sequence = []
    for part in new_state:
        new_positions.append((part[0], part[1]))
        new_sequence.append(part[2])

    reward1 = reward_function(new_positions, new_sequence)
    reward0 = reward_function(positions, sequence)
    
    delta = -(reward1 - reward0)
    return delta


def policy_gradient_main_loop(protein_sequence, num_episodes, learning_rate): ## select action ->heuristiek
    for episode in range(num_episodes):
        current_state = protein_sequence
        episode_data = []

        for step in range(100):
            action = select_action(current_state, policy_weights) # Alleen nog heuristic
            new_state = apply_action(current_state, action) ## Check
            reward = compute_reward(new_state, current_state) #  Check
            episode_data.append((current_state, action, reward))

            current_state = new_state

        policy_weights = update_policy_weights(episode_data, policy_weights, learning_rate)

    return policy_weights

##################ACTION#############################################
def select_action(state, policy_weights):
    possible_actions = get_possible_actions(state) #check
    action_scores = []

    for action in possible_actions:
        score = compute_heuristic_score(state, action, policy_weights)
        action_scores.appemd(score)

    max_score = max(action_scores)
    best_actions = [action for action, score in zip(possible_actions,action_scores) if score == max_score]
    selected_action = random.choice(best_actions)

    return selected_action

def apply_action(state, index, clockwise=True): ## check
    current_state = [list(sequence) for sequence in state]  # Convert tuples to lists for mutability
    positions = [(sequence[0], sequence[1]) for sequence in state]

    if index <= 0 or index >= len(positions) - 1:
        return current_state

    pivot = positions[index]
    new_positions = positions.copy()
    for i in range(index + 1, len(positions)):
        dx, dy = positions[i][0] - pivot[0], positions[i][1] - pivot[1]
        if clockwise:
            new_positions[i] = (pivot[0] - dy, pivot[1] + dx)
        else:
            new_positions[i] = (pivot[0] + dy, pivot[1] - dx)

    if is_valid_configuration(new_positions):   
        for i, position in enumerate(new_positions):
            current_state[i][0] = position[0]
            current_state[i][1] = position[1]
        return current_state

    return current_state

def validate_action(state, action, clockwise=True): ##check
    index = action[0]
    clockwise = action[1]
    current_state = [list(sequence) for sequence in state]  # Convert tuples to lists for mutability
    positions = [(sequence[0], sequence[1]) for sequence in state]

    if index <= 0 or index >= len(positions) - 1:
        return False

    pivot = positions[index]
    new_positions = positions.copy()
    for i in range(index + 1, len(positions)):
        dx, dy = positions[i][0] - pivot[0], positions[i][1] - pivot[1]
        if clockwise:
            new_positions[i] = (pivot[0] - dy, pivot[1] + dx)
        else:
            new_positions[i] = (pivot[0] + dy, pivot[1] - dx)

    if is_valid_configuration(new_positions):   
        for i, position in enumerate(new_positions):
            current_state[i][0] = position[0]
            current_state[i][1] = position[1]
        return True

    return False

def get_possible_actions(state): ##check
    list = []
    for i, part  in enumerate(state):
        actiontrue = (i, True)
        actionfalse = (i, False)
        if validate_action(state, actiontrue):
            list.append(actiontrue)
        if validate_action(state, actionfalse):
            list.append(actionfalse)
    return list


def is_valid_configuration(positions): ##check
    return len(positions) == len(set(positions))


##############################################################
def update_policy_weights(episode_data, current_weights, learning_rate):
    updated_weights = current_weights.copy()

    for state, action, reward in episode_data:
        for i in range(len(current_weights)):
            heuristic_influence = estimate_heuristic_influence(i, state, action)

            if reward > 0 and heuristic_influence > 0 or reward < 0 and heuristic_influence  < 0:
                updated_weights[i] += learning_rate 
            elif reward > 0 and heuristic_influence < 0 or reward < 0 and heuristic_influence > 0:
                updated_weights -= learning_rate

    return updated_weights

############################## heuristics #################
def compute_heuristic_score(state, action, policy_weights):
    total = 1
    return total

###########################################################