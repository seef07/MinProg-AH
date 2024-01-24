import random
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

sequence = "HHPHHHPHPHHHPH"

energy_matrix = {
    ##
}

policy_weights = [0.2,0.2,0.2,0.2,0.2]

episodes = 1000

data = [(i, 0, amino_acid) for i, amino_acid in enumerate(sequence)]

print(data[1])

def dompute_heuristic_score(state, action, policy_weights):
    total = 1
    return total

def reward_function(positions, sequence):
    energy = 0
    for i in range(len(sequence)):
        for j in range(i + 2, len(sequence)):  # Start from i+2 to skip consecutive amino acids
            if abs(positions[i][0] - positions[j][0]) + abs(positions[i][1] - positions[j][1]) == 1:
                pair = ''.join(sorted([sequence[i], sequence[j]]))
                energy += energy_matrix[pair]
    return energy

def compute_reward(new_state, current_state):
    ##parameter preperation
    ## return new - current
    return True

def policy_gradient_main_loop(protein_sequence, num_episodes, learning_rate):
    for episode in range(num_episodes):
        current_state = protein_sequence
        episode_data = []

        for step in range(100):
            action = select_action(current_state, policy_weights)
            new_state = apply_action(current_state, action)
            reward = compute_reward(new_state, current_state)
            episode_data.append((current_state, action, reward))

            curremt_state = new_state

        policy_weights = update_policy_weights(episode_data, policy_weights, learning_rate)

    return policy_weights

def select_action():
    possible_actions = get_possible_actions(state)
    action_scores = []

    for action in possible_actions:
        score = compute_heuristic_score(state, action, policy_weights)
        action_scores.appemd(score)

    max_score = max(action_scores)
    best_actions = [action in action, score in zip(possible_actions,action_scores) if score == max_score]
    selected_action = random.choice(best_actions)
    return select_action

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