import random
import pdb
import sys
sys.path.append('..')
from tqdm import tqdm
from typing import TYPE_CHECKING
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score
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
    def __init__(self, protein, dimensions, learning_rate=0.1, discount_factor=0.9, epsilon=1.0, epsilon_min=0.01, epsilon_decay=0.995):
        super().__init__(protein, dimensions)
        self.q_table = QTable()
        self.learning_rate = learning_rate   # Set learning_rate as an instance attribute
        self.discount_factor = discount_factor  # Set discount_factor as an instance attribute
        self.epsilon = epsilon
        self.epsilon_min = epsilon_min
        self.epsilon_decay = epsilon_decay


    def softmax_probabilities(self, q_values):
        exp_q = np.exp(q_values - np.max(q_values))  # Subtract max to prevent overflow
        return exp_q / exp_q.sum()

    def choose_action(self, state, possible_actions, max_attempts=10):
        attempts = 0
        while attempts < max_attempts:
            if np.random.rand() < self.epsilon:
                # Exploration: choose a random action
                action = random.choice(possible_actions)
            else:
                # Exploitation: choose the best action based on Q-values
                q_values = [self.q_table.get(state, action) for action in possible_actions]
                max_q_value = max(q_values)
                best_actions = [action for action, q in zip(possible_actions, q_values) if q == max_q_value]
                action = random.choice(best_actions)

            # Simulate new order and check validity
            new_order = self.simulate_new_order(state, action)
            if fast_validate_protein(new_order):
                return action

            attempts += 1

        # Fallback strategy if no valid action is found
        return random.choice(possible_actions)  
    def simulate_new_order(self, state, action):
        # Create a copy of the current order
        new_order = list(state.get_order())

        # Update the new order based on the chosen action
        node_idx, new_direction = action
        new_order[node_idx - 1] = new_direction  # Adjust indexing if necessary

        return new_order
    
    def learn(self, state, action, reward, next_state, possible_actions):
        self.q_table.update(state, action, reward, self.learning_rate, self.discount_factor, next_state, possible_actions)
        self.epsilon = max(self.epsilon * self.epsilon_decay, self.epsilon_min)

def get_free_directions(current_state, directions):
    free_actions = []
    for idx, node in enumerate(current_state.nodes[1:], start=1):  # start=1 to account for skipping first node
        node_free_directions = node.get_free_directions(directions)
        node_actions = [(idx, direction) for direction in node_free_directions]
        free_actions.extend(node_actions)
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

def get_reward(current_state, next_state):
    max_score = 40
        # Define the weights for different components of the reward
    score_weight =1.0  # Adjust as necessary
    fold_amount_weight = 1  # Adjust as necessary
    dimension_penalty_weight = 1  # Adjust as necessary

    if not fast_validate_protein(next_state.get_order()):
        return -1  # Penalty for invalid states

    current_score = -current_state.get_bond_score() if fast_validate_protein(current_state.get_order()) else 0
    next_score = -next_state.get_bond_score()

    # Normalize the scores
    normalized_current_score = current_score / max_score
    normalized_next_score = next_score / max_score

    # Calculate fold amount
    current_fold_amount = calculate_fold_amount(current_state)
    next_fold_amount = calculate_fold_amount(next_state)
    fold_amount_change = next_fold_amount - current_fold_amount

    # Dimension penalty
    dimension_penalty = get_dimension_penalty(next_state)

    # Compute the final reward
    reward = (normalized_next_score - normalized_current_score) * score_weight
    reward += fold_amount_change * fold_amount_weight
    reward -= dimension_penalty * dimension_penalty_weight

    return reward


def run_protein_folding(sequence, iteration):
    protein = Protein(sequence)
    agent = ProteinFoldingAgent(protein, dimensions = 2, learning_rate=0.1, discount_factor=0.9, epsilon=0.9, epsilon_min=0.01, epsilon_decay=0.995)
    num_iterations = iteration

    current_state = protein

    for _ in range(num_iterations):
        possible_actions = get_free_directions(current_state, agent.directions)
        chosen_action = agent.choose_action(current_state, possible_actions)
        next_state = get_next_state(current_state, chosen_action)
        reward = get_reward(current_state, next_state)

        agent.learn(current_state, chosen_action, reward, next_state, possible_actions)

        print(f'Current State: {current_state.get_order()}, Action: {chosen_action}, Next State: {next_state}, Reward: {reward}/{current_state.get_bond_score()}')

        current_state = next_state

        if current_state == 'goal_state':
            break

    return protein
