import random
import pdb
import sys
sys.path.append('..')
from tqdm import tqdm
from typing import TYPE_CHECKING
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score, fast_validate_protein_action
from protein_folding.protein import Protein

from .. import Algorithm
from ..heuristics import *
import numpy as np

class QTable:
    def __init__(self):
        self.table = {}

    def get(self, state, action):
        return self.table.get((state, action), 0)

    def set(self, state, action, value):
        self.table[(state, action)] = value

    def update(self, state, action, reward, learning_rate, discount_factor, next_state, possible_actions):
        old_value = self.get(state, action)
        future_rewards = [self.get(next_state, next_action) for next_action in possible_actions]
        max_future_reward = max(future_rewards, default=0)
        new_value = old_value + learning_rate * (reward + discount_factor * max_future_reward - old_value)
        self.set(state, action, new_value)

class ProteinFoldingAgent(Algorithm):
    def __init__(self, protein, dimensions, learning_rate=0.1, discount_factor=0.9, epsilon=0.5, epsilon_min=0.01, epsilon_decay=0.995):
        super().__init__(protein, dimensions)
        self.q_table = QTable()
        self.learning_rate = learning_rate   # Set learning_rate as an instance attribute
        self.discount_factor = discount_factor  # Set discount_factor as an instance attribute
        self.epsilon = epsilon
        self.epsilon_min = epsilon_min
        self.epsilon_decay = epsilon_decay
        self.protein_sequence = protein.sequence
    	

    def softmax_probabilities(self, q_values):
        exp_q = np.exp(q_values - np.max(q_values))  # Subtract max to prevent overflow
        return exp_q / exp_q.sum()

    def choose_action(self, state, possible_actions):
        if not possible_actions:
            raise Exception("No valid actions available")
        if random.random() < self.epsilon:
            # Exploration: choose a random action
            chosen_action = random.choice(possible_actions)
        else:
            # Exploitation: choose the best action based on Q-table values
            q_values = [self.q_table.get(state, action) for action in possible_actions]
            max_q_value = max(q_values)
            # In case there are multiple actions with the same max Q-value, choose randomly among them
            best_actions = [action for action, q_value in zip(possible_actions, q_values) if q_value == max_q_value]
            chosen_action = random.choice(best_actions)

        # Decay the epsilon value to gradually shift from exploration to exploitation
        self.epsilon = max(self.epsilon * self.epsilon_decay, self.epsilon_min)
        print(chosen_action)
        return chosen_action

    def simulate_new_order(self, state, action):
        # Create a copy of the current order
        new_order = list(state.get_order())

        # Update the new order based on the chosen action
        node_idx, new_direction = action
        new_order[node_idx] = new_direction  # Adjust indexing if necessary

        return new_order
    
    def learn(self, state, action, reward, next_state, possible_actions):
        self.q_table.update(state, action, reward, self.learning_rate, self.discount_factor, next_state, possible_actions)
        self.epsilon = max(self.epsilon * self.epsilon_decay, self.epsilon_min)

def get_free_directions(current_state, directions):
    current = current_state.get_order()
    free_actions = []
    for idx in range(2, len(current)):  # start=1 to account for skipping first node
        for direction in directions:
            vector = (idx, direction)
            if fast_validate_protein_action(current, vector):
                free_actions.append(vector)      
    return free_actions  # Return a list of (node_idx, direction) tuples

def get_next_state(current_state, action):
    # Logic to determine the next state based on the current state and chosen action
    node_idx, new_direction = action  # Action decomposed into node index and new direction
    current_state.nodes[node_idx].change_direction(new_direction)  
    return current_state

def calculate_fold_amount(state):
    # Create an instance of FoldAmount with the given state
    fold_amount_heuristic = FoldAmount(state)
    
    # Run the heuristic to calculate scores
    score = fold_amount_heuristic.calculate(state)

    return score

def are_states_same(current_state, next_state):
    # Assuming the state comparison is based on the order of nodes
    return current_state.get_order() == next_state.get_order()

def get_dimension_penalty(state):
    dimension_heuristic = MinimiseDimensions(state)
    dimension_heuristic.run()
    scores = dimension_heuristic.interpret()
    return sum(scores)

def get_reward(current_state, next_state, action, possible_actions, directions):
    reward = 0
    current_score = -current_state.get_bond_score() if fast_validate_protein(current_state.get_order()) else 0
    next_score = -next_state.get_bond_score()
    #total_reward = lookahead_reward(current_state, action, directions, max_depth = 3, current_depth = 0)
    next_fold_amount = calculate_fold_amount(next_state)
    # Compute the final reward
    dimension_penalty = get_dimension_penalty(next_state)*10
    print(f"{reward}+{next_score}++{next_fold_amount}")
    reward += next_score*2 + next_fold_amount /4  
    return reward

def lookahead_reward(state, action, directions, max_depth, current_depth=0):
    if current_depth == max_depth:
        return 0  # No further rewards at maximum depth

    next_state = get_next_state(state, action)
    immediate_reward = -next_state.get_bond_score()
    lookahead_rewards = 0
    print(state.get_order())
    possible_actions = get_free_directions(next_state, directions)
    
    for next_action in possible_actions:
        print(f"Nextaction{next_action}")
        lookahead_rewards += lookahead_reward(next_state, next_action, directions, max_depth, current_depth + 1)
    
    return immediate_reward + lookahead_rewards


def run_protein_folding(sequence, iteration):
    protein = Protein(sequence)
    agent = ProteinFoldingAgent(protein, dimensions = 2, learning_rate=0.4, discount_factor=0.9, epsilon=0.3, epsilon_min=0.01, epsilon_decay=0.995)
    num_iterations = iteration

    current_state = protein
    best_state = 0
    for _ in range(num_iterations):
        possible_actions = get_free_directions(current_state, agent.directions)
        print(possible_actions)
        chosen_action = agent.choose_action(current_state, possible_actions)
        next_state = get_next_state(current_state, chosen_action)
        reward = get_reward(current_state, next_state, chosen_action, possible_actions,  agent.directions)

        agent.learn(current_state, chosen_action, reward, next_state, possible_actions)

        getbond = current_state.get_bond_score()
        best_score = 0
        if getbond <= best_score and fast_validate_protein(current_state.get_order()):
            best_state = current_state
            best_score = getbond
            print("#####################################!!!################################")

        print(f'Current State: {current_state.get_order()}, Action: {chosen_action}, Next State: {next_state}, Reward: {reward}/{current_state.get_bond_score()}')

        current_state = next_state

        if current_state == 'goal_state':
            break


    return best_state, current_state
