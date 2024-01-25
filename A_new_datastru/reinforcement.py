import random
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from heuristics.heuristics import altheuristic, currentletter, compactness_heuristic, folding_heuristic

sequence = "HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP"

energy_matrix = {
    'HH': -1, 'CC': -5, 'CH': -1, 'HC':-1, 
    'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
}


policy_weights = [-0.27 , -0.122, -0.035, -0.07, 0.13]



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

scores_per_iteration = []
def policy_gradient_main_loop(protein_sequence, num_episodes, learning_rate, policy_weights):
    global scores_per_iteration 
    bestscore = 0
    beststate = protein_sequence 
    i = 0
    for episode in range(num_episodes):
        current_state = protein_sequence
        episode_data = []
        print("start")
        for step in range(10000):
            i += 1
            action = select_action(current_state, policy_weights)  # Use heuristic
            new_state = apply_action(current_state, action[0], action[1])  # Apply action
            reward = compute_reward(new_state, current_state)  # Compute reward
            bondscore = reward_function(extract_positions(new_state), sequence)
            if bestscore > bondscore:
                bestscore = bondscore
                beststate = new_state
            episode_data.append((current_state, action, reward))
            scores_per_iteration.append((i, bondscore))
            current_state = new_state

        print("#################################################################################################")
        print(f'bestscore: {bestscore} Beststate: {beststate}')
        policy_weights = update_policy_weights(episode_data, policy_weights, learning_rate)
        print(policy_weights)
    return policy_weights

def extract_positions(state):
    return [(x, y) for x, y, _ in state]

iteration = 0
##################ACTION#############################################
def select_action(state, policy_weights, start_temp=1.0, end_temp=0.01):
    global iteration
    iteration += 1
    temp = start_temp - iteration * (start_temp - end_temp) / 100000
    possible_actions = get_possible_actions(state) #check
    action_scores = []

    for action in possible_actions:
        score = compute_heuristic_score(state, action, policy_weights)
        action_scores.append(score)

    max_score = max(action_scores)
    best_actions = [action for action, score in zip(possible_actions,action_scores) if score == max_score]
    selected_action = random.choice(best_actions)
    new_state = apply_action(state, selected_action[0], selected_action[1])  
    new_energy = reward_function(extract_positions(new_state), sequence)
    current_energy = reward_function(extract_positions(state), sequence)
    exp_argument = (new_energy - current_energy) / temp
    # Cap the argument to avoid overflow
    capped_argument = min(exp_argument, 1)  # 20 is an arbitrary cap, adjust as needed

    if random.random() < math.exp(capped_argument):
        selected_action = random.choice(possible_actions)

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
########################


##############################################################
def update_policy_weights(episode_data, current_weights, learning_rate):
    updated_weights = current_weights.copy()

    for state, action, reward in episode_data:
        for i in range(len(current_weights)):
            heuristic_influence = estimate_heuristic_influence(i, state, action)

            if reward > 0 and heuristic_influence > 0 or reward < 0 and heuristic_influence  < 0:
                updated_weights[i] += learning_rate 
            elif reward > 0 and heuristic_influence < 0 or reward < 0 and heuristic_influence > 0:
                updated_weights[i] -= learning_rate
    return updated_weights

############################## heuristics #################
def compute_heuristic_score(state, action, policy_weights):
    nextstate = apply_action(state, action[0], action[1])

    altscore = altheuristic(sequence, action[0])
    ptscore = currentletter(sequence, action[0])
    compactscore = compactness_heuristic(state, nextstate)
    patternscore = folding_heuristic(sequence, action[0])
    rewardscore = compute_reward(nextstate, state)

    # Debugging: Print the heuristic scores

    heuristic_scores = [altscore, ptscore, compactscore, patternscore, rewardscore]

    total = sum(weight * score for weight, score in zip(policy_weights, heuristic_scores))

    return total


###########################################################

def estimate_heuristic_influence(i, state, action):
    nextstate = apply_action(state, action[0], action[1])

    if i == 0:
        return altheuristic(sequence, action[0])
    if i == 1:
        return currentletter(sequence, action[0])
    if i == 2:
        return compactness_heuristic(state, nextstate)
    if i == 3:
        return folding_heuristic(sequence, action[0])
    if i == 4:
        return compute_reward(nextstate, state)
    
import matplotlib.pyplot as plt

def visualize_protein(state):
    # Extract positions and sequence from the state
    positions = [(x, y) for x, y, _ in state]
    sequence = [acid for _, _, acid in state]

    # Map each amino acid type to a color
    colors = {'H': 'red', 'P': 'blue', 'C': 'green'}
    color_sequence = [colors[acid] for acid in sequence]

    x, y = zip(*positions)
    plt.figure(figsize=(10, 10))
    plt.scatter(x, y, color=color_sequence)

    for i in range(len(positions) - 1):
        plt.plot([positions[i][0], positions[i + 1][0]], [positions[i][1], positions[i + 1][1]], color='black')

    plt.title(f"Protein Folding Visualization\nEnergy: {reward_function(positions, sequence)}")
    plt.show()

# Example usage
state =   [[0, 0, 'H'], [1, 0, 'H'], [1, 1, 'P'], [2, 1, 'C'], [3, 1, 'H'], [4, 1, 'H'], [4, 0, 'P'], [3, 0, 'C'], [2, 0, 'C'], [2, -1, 'P'], [3, -1, 'C'], [3, -2, 'P'], [2, -2, 'P'], [1, -2, 'H'], [1, -1, 'H'], [0, -1, 'H'], [0, -2, 'H'], [0, -3, 'P'], [-1, -3, 'P'], [-1, -2, 'H'], [-1, -1, 'C'], [-1, 0, 'H'], [-1, 1, 'P'], [-1, 2, 'H'], [0, 2, 'P'], [1, 2, 'H'], [2, 2, 'C'], [3, 2, 'H'], [3, 3, 'P'], [4, 3, 'P']]
visualize_protein(state)



print(policy_gradient_main_loop(protein_sequence = data, num_episodes=15, learning_rate= 0.001, policy_weights= policy_weights))

steps, bondscores = zip(*scores_per_iteration)

# Create a plot
plt.figure(figsize=(10, 6))
plt.plot(steps, bondscores, marker='o')  # 'o' creates circular markers on each data point
plt.title('Bond Scores Per Iteration')
plt.xlabel('Iteration Step')
plt.ylabel('Bond Score')
plt.grid(True)
plt.show()