
import random
import math
from heuristics.heuristics import altheuristic, currentletter, compactness_heuristic, folding_heuristic, distance_heuristic, compactness_heuristic_H, compactness_heuristic_C
import matplotlib.pyplot as plt
import numpy as np 

EPISODE = 0

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
        self.weights_history = []  # Store weights history
        self.bond_scores_history = [] 

    def _initialize_state(self):
        return [(i, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]

    def run(self, num_episodes, learning_rate):
        print("Start run")
        global iteration
        global EPISODE
       
        randomize_weights_interval = 30
        for episode in range(num_episodes):
            EPISODE = 300
            print(episode)
            iteration = 0
            self._run_episode(learning_rate)

            if randomize_weights_interval and episode % randomize_weights_interval == 0 and episode != 0:
                self.randomize_weights()


    def _run_episode(self, learning_rate):
        global EPISODE
        current_state = [(i, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]
        state_history = set()
        last_action = None
        action_repetition_count = 0
        repetition_threshold = 45  
        actions2 = None
        episode_data = []

        for step in range(EPISODE):
            action, heuristic_influence = self.action_selector.select_action(current_state)
            state_instance = State(current_state)
            new_state = state_instance.apply_action(action[0], action[1])  # Apply action
            #self.visualize_in_cmd(new_state)
            reward = state_instance.compute_reward(new_state)  # Compute reward
            bondscore = state_instance.reward_function()
            self.bond_scores_history.append(bondscore)
            new_state_instance = State(new_state)
            new_bondscore = new_state_instance.reward_function()
            self.weights_history.append(self.policy_weights.copy())
            if self.bestscore > new_bondscore:
                self.bestscore = new_bondscore
                self.beststate = new_state
            
            episode_data.append((current_state, action, reward, heuristic_influence))
            current_state = new_state

            if action == last_action or action == actions2:
                action_repetition_count += 1
                if action_repetition_count > repetition_threshold:
                    print("###########################################reset#######")
                    self.randomize_weights()
                    current_state = [(i, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]
                    action_repetition_count = 0
            else:
                action_repetition_count = 0
            actions2 = last_action
            last_action = action
            


        print("###############################")
        print(self.bestscore)
        print(self.policy_weights)
        print(self.beststate)
        print("###############################")
        self.policy_weights = self.update_policy_weights(self.policy_weights, episode_data, learning_rate)



    def randomize_weights(self):
        # Randomize the policy weights by adding a random value between -5 and 5
        self.policy_weights = [np.random.uniform(0.01, 1) for w in self.policy_weights]

    @staticmethod
    def update_policy_weights(current_weights, episode_data, learning_rate):
        updated_weights = current_weights.copy()

        for state, action, reward, heuristic_influence in episode_data:
            for i in range(len(current_weights)):
                # Adjust weights based on reward and heuristic influence
                if reward > 0 and heuristic_influence[i] > 0 or reward < 0 and heuristic_influence[i] < 0:
                    updated_weights[i] *= (1+learning_rate/2) 
                elif reward > 0 and heuristic_influence[i] < 0 or reward < 0 and heuristic_influence[i] > 0:
                    updated_weights[i] *= (1-learning_rate/3)
            #    elif reward == 0 and heuristic_influence[i] > 0 or reward == 0 and heuristic_influence[i] > 0:
            #       updated_weights[i] *= (1-learning_rate/4)
                
                # Clip the updated weight to ensure it stays within -1 and 1
                updated_weights[i] = np.clip(updated_weights[i], 0, 1)
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

    
    def plot_weights_and_bond_scores(self):
        # Plotting weights
        weights_array = np.array(self.weights_history)
        plt.figure(figsize=(12, 6))
        for i in range(weights_array.shape[1]):
            plt.plot(weights_array[:, i], label=f'Weight {i+1}')
        plt.plot(self.bond_scores_history, label='Bond Score')
        plt.title('Policy Weights Over Iterations')
        plt.xlabel('Iteration')
        plt.ylabel('Weight Value')
        plt.legend()
        plt.show()

    def visualize_in_cmd(self, new_state):
        positions = [(x, y) for x, y, _ in new_state]
        sequence = [acid for _, _, acid in new_state]
        # Determine the size of the grid
        x_coords, y_coords = zip(*positions)
        min_x, max_x = min(x_coords), max(x_coords)
        min_y, max_y = min(y_coords), max(y_coords)

        # Adjust grid size for connections
        grid_width = (max_x - min_x + 1) * 2 + 1
        grid_height = (max_y - min_y + 1) * 2 + 1

        # Create an empty grid
        grid = [[' ' for _ in range(grid_width)] for _ in range(grid_height)]

        # Place the amino acids and connections on the grid
        for i, (pos, acid) in enumerate(zip(positions, sequence)):
            x, y = (pos[0] - min_x) * 2 + 1, (pos[1] - min_y) * 2 + 1
            grid[y][x] = acid

            # Draw connections
            if i > 0:
                prev_x, prev_y = (positions[i-1][0] - min_x) * 2 + 1, (positions[i-1][1] - min_y) * 2 + 1
                if prev_x == x:  # Vertical connection
                    grid[min(y, prev_y) + 1][x] = '|'
                elif prev_y == y:  # Horizontal connection
                    grid[y][min(x, prev_x) + 1] = '-'

        # Print the grid
        for row in grid:
            print(''.join(row))


iteration = 0
class ActionSelector:
    def __init__(self, policy_weights):
        self.policy_weights = policy_weights

    def select_action(self, state, start_temp=0.1, end_temp=0.01):
        global EPISODE 
        global iteration
        iteration += 1
        temp = start_temp - iteration * (start_temp - end_temp) / EPISODE
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


        if random.random() < 0.1:
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
        distance = distance_heuristic(state, nextstate)
        hcompact = compactness_heuristic_H(state, nextstate)
        ccompac = compactness_heuristic_C(state, nextstate)

        heuristic_scores = [ccompac, ptscore, compactscore, hcompact, rewardscore, distance]

        # Calculate the total score using policy weights
        total = sum(weight * score for weight, score in zip(policy_weights, heuristic_scores))
        new_state = State(nextstate)
        total -= new_state.reward_function()
        # Create a vector indicating the sign of each heuristic score
        sign_vector = [1 if score > 0 else -1 if score < 0 else 0 for score in heuristic_scores]

        return total, sign_vector




#update policy weight
sequence = "HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH"
policy_weights = [0.2, 0.2, 0.2, 0.2, 0.2, 0.1]

simulator = ProteinFoldingSimulator(sequence, policy_weights)
simulator.run(num_episodes=500, learning_rate=0.02)

print(simulator.bestscore)
print(simulator.beststate)
simulator.visualize_protein()
simulator.plot_weights_and_bond_scores()