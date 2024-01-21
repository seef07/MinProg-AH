import random
import pdb
import sys
sys.path.append('..')
from tqdm import tqdm
from typing import TYPE_CHECKING
from protein_folding.fast_protein import fast_validate_protein, fast_compute_bond_score
from protein_folding.protein import Protein

from .. import Algorithm


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
        self.epsilon = epsilon
        self.epsilon_min = epsilon_min
        self.epsilon_decay = epsilon_decay

    def choose_action(self, state, possible_actions):
        if random.random() < self.epsilon:
            return random.choice(possible_actions)
        else:
            q_values = [self.q_table.get(state, action) for action in possible_actions]
            max_q_value = max(q_values)
            max_actions = [action for action, q in zip(possible_actions, q_values) if q == max_q_value]
            return random.choice(max_actions) if max_actions else random.choice(possible_actions)

    def learn(self, state, action, reward, next_state, possible_actions):
        self.q_table.update(state, action, reward, self.learning_rate, self.discount_factor, next_state, possible_actions)
        self.epsilon = max(self.epsilon * self.epsilon_decay, self.epsilon_min)

def get_free_directions(current_state, directions):
    # and each node in this structure can provide its free directions
    free_directions = []
    for node in current_state.nodes[1:]:
        node_free_directions = node.get_free_directions(directions)  
        free_directions.extend(node_free_directions)
    return list(set(free_directions))  # Return unique free directions

def get_next_state(current_state, action):
    # Logic to determine the next state based on the current state and chosen action
    node_idx, new_direction = action  # Action decomposed into node index and new direction
    new_state = current_state.set_new_direction(node_idx, new_direction)  # Method to update the state
    return new_state

def get_reward(current_state, next_state):
    # Logic to calculate the reward based on the current and next states
    # Assuming reward is based on the bond score
    if fast_validate_protein(current_state):
        current_score = current_state.get_bond_score()  # Method to get bond score of the current state
    if fast_validate_protein(next_state):
        next_score = next_state.get_bond_score()  # Method to get bond score of the next state
    reward = next_score - current_score  # Reward is the improvement in bond score
    return reward

def run_protein_folding(sequence, iteration):
    protein = Protein(sequence)
    agent = ProteinFoldingAgent(protein, dimensions = 2, learning_rate=0.1, discount_factor=0.9, epsilon=1.0, epsilon_min=0.01, epsilon_decay=0.995)
    num_iterations = iteration

    current_state = protein

    for _ in range(num_iterations):
        possible_actions = get_free_directions(current_state, agent.directions)
        chosen_action = agent.choose_action(current_state, possible_actions)
        next_state = get_next_state(current_state, chosen_action)
        reward = get_reward(current_state, next_state)

        agent.learn(current_state, chosen_action, reward, next_state, possible_actions)

        print(f'Current State: {current_state}, Action: {chosen_action}, Next State: {next_state}, Reward: {reward}')

        current_state = next_state

        if current_state == 'goal_state':
            break
