import random
import math
from heuristics.heuristics import altheuristic, currentletter, compactness_heuristic, folding_heuristic
import matplotlib.pyplot as plt
import numpy as np 
class ProteinFoldingSimulator:
    def __init__(self, sequence, policy_weights):
        self.sequence = sequence
        self.energy_matrix = {
                                'HH': -1, 'CC': -5, 'CH': -1, 'HC':-1, 
                                'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
                            }
        self.policy_weights = policy_weights
        self.action_selector = ActionSelector(policy_weights)
        self.state = self._initialize_state()
        self.scores_per_iteration = []
        self.beststate = None
        self.bestscore = 0

    def _initialize_state(self):
        return [(i, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]

    def run(self, num_episodes, learning_rate):
        randomize_weights_interval = 5
        for episode in range(num_episodes):
            self._run_episode(learning_rate)

            if randomize_weights_interval and episode % randomize_weights_interval == 0 and episode != 0:
                self.randomize_weights()


    def _run_episode(self, learning_rate):
        current_state = self.state
        episode_data = []
        for step in range(2000):
            action, heuristic_influence = self.action_selector.select_action(current_state)
            state_instance = State(current_state)
            new_state = state_instance.apply_action(action[0], action[1])  # Apply action
            reward = state_instance.compute_reward(new_state)  # Compute reward
            bondscore = state_instance.reward_function()
            if self.bestscore > bondscore:
                self.bestscore = bondscore
                self.beststate = new_state
            
            episode_data.append((current_state, action, reward, heuristic_influence))
            current_state = new_state
        print(self.bestscore)
        self.policy_weights = self.update_policy_weights(self.policy_weights, episode_data, learning_rate)

    def randomize_weights(self, scale=0.1):
        # Randomize the policy weights by adding a small random value and then clipping
        self.policy_weights = [np.clip(w + np.random.uniform(-scale, scale), -1, 1) for w in self.policy_weights]

    @staticmethod
    def update_policy_weights(current_weights, episode_data, learning_rate):
        updated_weights = current_weights.copy()

        for state, action, reward, heuristic_influence in episode_data:
            for i in range(len(current_weights)):
                # Adjust weights based on reward and heuristic influence
                if reward > 0 and heuristic_influence[i] > 0 or reward < 0 and heuristic_influence[i] < 0:
                    updated_weights[i] += learning_rate 
                elif reward > 0 and heuristic_influence[i] < 0 or reward < 0 and heuristic_influence[i] > 0:
                    updated_weights[i] -= learning_rate
                
                # Clip the updated weight to ensure it stays within -1 and 1
                updated_weights[i] = np.clip(updated_weights[i], -1, 1)

        return updated_weights


    def visualize_protein(self):
        state = self.beststate
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

        plt.title(f"Protein Folding Visualization\nEnergy: {self.bestscore}")
        plt.show()

iteration = 0
class ActionSelector:
    def __init__(self, policy_weights):
        self.policy_weights = policy_weights

    def select_action(self, state, start_temp=1.0, end_temp=0.01):
        global iteration
        iteration += 1
        temp = start_temp - iteration * (start_temp - end_temp) / 2000
        possible_actions = self._get_possible_actions(state)
        action_scores = []

        for action in possible_actions:
            # Gebruik self.policy_weights
            score, signvector = Heuristic.compute_heuristic_score(state, action, self.policy_weights) 
            action_scores.append(score)

        max_score = max(action_scores)
        best_actions = [action for action, score in zip(possible_actions,action_scores) if score == max_score]
        selected_action = random.choice(best_actions)
        state_instance = State(state)
        new_state = state_instance.apply_action(selected_action[0], selected_action[1])
        new_state_instance = State(new_state)
        new_energy = new_state_instance.reward_function()
        current_energy = state_instance.reward_function()     
        exp_argument = (new_energy - current_energy) / temp

        capped_argument = min(exp_argument, 1) 

        if random.random() < math.exp(capped_argument):
            selected_action = random.choice(possible_actions)

        score, heuristic_influence = Heuristic.compute_heuristic_score(state, selected_action, policy_weights) 
        return selected_action, heuristic_influence

    def _get_possible_actions(self, state): 
        list = []
        for i, part  in enumerate(state):
            actiontrue = (i, True)
            actionfalse = (i, False)
            if self._validate_action(state, actiontrue):
                list.append(actiontrue)
            if self._validate_action(state, actionfalse):
                list.append(actionfalse)
        return list
    
    def _validate_action(sself, state, action, clockwise=True): 
        state_class = State(state)
        index = action[0]
        clockwise = action[1]
        current_state = [list(sequence) for sequence in state]
        positions = state_class.positions

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
        if state_class.is_valid_configuration(new_positions):   
            for i, position in enumerate(new_positions):
                current_state[i][0] = position[0]
                current_state[i][1] = position[1]
            return True

        return False


class State:
    def __init__(self, state):
        self.positions = [(x, y) for x, y, _ in state]
        self.sequence = [z for _, _, z in state]
        self.state = state
        self.energy_matrix = {
                                'HH': -1, 'CC': -5, 'CH': -1, 'HC':-1, 
                                'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
                            }
    def apply_action(self, index, clockwise=True):
        current_state = [list(sequence) for sequence in self.state]  # Convert tuples to lists for mutability
        positions = self.positions

        # Check if the action is valid; if not, return the original state immediately
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

        # Check if the new configuration is valid
        if not self.is_valid_configuration(new_positions):
            return current_state

        # Update the state with new positions
        for i, position in enumerate(new_positions):
            current_state[i][0] = position[0]
            current_state[i][1] = position[1]
        return current_state
    
    
    def get_positions(self):
        return [(x, y) for x, y, _ in self.state]
    
    @staticmethod
    def is_valid_configuration(positions): 
        return len(positions) == len(set(positions))
    
    def reward_function(self):
        energy = 0
        for i in range(len(self.sequence)):
            for j in range(i + 2, len(self.sequence)):  # Start from i+2 to skip consecutive amino acids
                if abs(self.positions[i][0] - self.positions[j][0]) + abs(self.positions[i][1] - self.positions[j][1]) == 1:
                    pair = ''.join(sorted([self.sequence[i], self.sequence[j]]))
                    energy += self.energy_matrix[pair]
        return energy
        

    def compute_reward(self, new_state):
        # Save the current state
        current_positions, current_sequence = self.positions, self.sequence

        # Update state to the new state
        self.positions = [(x, y) for x, y, _ in new_state]
        self.sequence = [z for _, _, z in new_state]

        # Compute the new reward
        new_reward = self.reward_function()

        # Restore the original state
        self.positions, self.sequence = current_positions, current_sequence

        # Compute the original reward
        original_reward = self.reward_function()

        # Compute the change in reward
        delta = new_reward - original_reward
        return delta

class Heuristic:
    def evaluate(self, state):
        pass

    def compute_heuristic_score(state, action, policy_weights):
        cstate = State(state)
        nextstate = cstate.apply_action(action[0], action[1])
        sequence = cstate.sequence

        # Compute individual heuristic scores
        altscore = altheuristic(sequence, action[0])
        ptscore = currentletter(sequence, action[0])
        compactscore = compactness_heuristic(state, nextstate)
        patternscore = folding_heuristic(sequence, action[0])
        rewardscore = cstate.compute_reward(nextstate)

        heuristic_scores = [altscore, ptscore, compactscore, patternscore, rewardscore]

        # Calculate the total score using policy weights
        total = sum(weight * score for weight, score in zip(policy_weights, heuristic_scores))

        # Create a vector indicating the sign of each heuristic score
        sign_vector = [1 if score > 0 else -1 if score < 0 else 0 for score in heuristic_scores]

        return total, sign_vector




##policy_gradient_main_loop

#select_action  ---> Heuristis etcetera



#update policy weight
sequence = "HHPCHHPCCPCPPHHHHPPHCHPHPHCHPP"
policy_weights = [0, 0, 0, 0, 0]

simulator = ProteinFoldingSimulator(sequence, policy_weights)
simulator.run(num_episodes=15, learning_rate=0.001)

print(simulator.bestscore)
print(simulator.beststate)
simulator.visualize_protein()
