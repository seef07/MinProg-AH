import random
import math
from heuristics.heuristics import altheuristic, currentletter, compactness_heuristic, folding_heuristic, distance_heuristic
import matplotlib.pyplot as plt
import numpy as np 

EPISODE = 1

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
        return [(i, 0, 0,  amino_acid) for i, amino_acid in enumerate(self.sequence)]

    def run(self, num_episodes, learning_rate):
        print("Start run")
        global iteration
        global EPISODE
        randomize_weights_interval = 20
        for episode in range(num_episodes):
            EPISODE = episode
            iteration = 0
            self._run_episode(learning_rate, episode)

            #if randomize_weights_interval and episode % randomize_weights_interval == 0 and episode != 0:
             #   self.randomize_weights()


    def _run_episode(self, learning_rate, episode):
        current_state = [(i, 0, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]
        episode_data = []
        for step in range(episode):
            action, heuristic_influence = self.action_selector.select_action(current_state)
            state_instance = State(current_state)
            new_state = state_instance.apply_action(action[0], action[1])  # Apply action
            reward = state_instance.compute_reward(new_state)  # Compute reward
            bondscore = state_instance.reward_function()
            self.bond_scores_history.append(bondscore)
            print(bondscore)
            self.weights_history.append(self.policy_weights.copy())
            if self.bestscore > bondscore:
                self.bestscore = bondscore
                self.beststate = new_state
            
            episode_data.append((current_state, action, reward, heuristic_influence))
            current_state = new_state
        print("###############################")
        print(self.bestscore)
        print(self.policy_weights)
        print(self.beststate)
        print("###############################")
        self.policy_weights = self.update_policy_weights(self.policy_weights, episode_data, learning_rate)

    def randomize_weights(self):
        # Randomize the policy weights by adding a random value between -5 and 5
        self.policy_weights = [w + np.random.uniform(-1, 1) for w in self.policy_weights]

    @staticmethod
    def update_policy_weights(current_weights, episode_data, learning_rate):
        updated_weights = current_weights.copy()

        for state, action, reward, heuristic_influence in episode_data:
            for i in range(len(current_weights)):
                # Adjust weights based on reward and heuristic influence
                if reward > 0 and heuristic_influence > 0 or reward < 0 and heuristic_influence < 0:
                    updated_weights[i] *= (1+learning_rate/2) 
                elif reward > 0 and heuristic_influence < 0 or reward < 0 and heuristic_influence > 0:
                    updated_weights[i] *= (1-learning_rate/3)
                # Clip the updated weight to ensure it stays within -1 and 1
                updated_weights[i] = np.clip(updated_weights[i], 0.1, 1)
        return updated_weights


    def visualize_protein(self):
        state = self.beststate

        # Extract positions and sequence from the state
        positions = [(x, y, z) for x, y, z, _ in state]
        sequence = [acid for _, _, _, acid in state]

        # Map each amino acid type to a color
        colors = {'H': 'red', 'P': 'blue', 'C': 'green'}
        color_sequence = [colors[acid] for acid in sequence]

        x, y, z = zip(*positions)

        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(x, y, z, c=color_sequence, marker='o')

        for i in range(len(positions) - 1):
            ax.plot([positions[i][0], positions[i + 1][0]], 
                    [positions[i][1], positions[i + 1][1]], 
                    [positions[i][2], positions[i + 1][2]], 
                    color='black')

        ax.set_title(f"Protein Folding Visualization\nEnergy: {self.bestscore}")
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

iteration = 0
class ActionSelector:
    def __init__(self, policy_weights):
        self.policy_weights = policy_weights

    def select_action(self, state, start_temp=0.25, end_temp=0.01):
        global iteration
        global EPISODE
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


        if random.random() < temp:
            selected_action = random.choice(possible_actions)

        score, heuristic_influence = Heuristic.compute_heuristic_score(state, selected_action, policy_weights) 
        return selected_action, heuristic_influence

    def _get_possible_actions(self, state):
        actions = []
        for i in range(len(state)):
            for rotation in ['clockwise_xy', 'counterclockwise_xy', 'clockwise_xz', 'counterclockwise_xz', 'clockwise_yz', 'counterclockwise_yz']:
                if self._validate_action(state, i, rotation):
                    actions.append((i, rotation))
        return actions


    def _validate_action(self, state, index, rotation):
        cstate = State(state)
        positions =  [(sequence[0], sequence[1], sequence[2]) for sequence in state]

        if index < 0 or index >= len(positions) - 1:
            return False

        pivot = positions[index]
        new_positions = positions.copy()

        dx, dy, dz = 0, 0, 0

        if rotation == 'clockwise_xy':
            # Clockwise rotation in the XY plane (horizontal plane)
            # Adjust the coordinates accordingly for a 3D rotation
            dx = pivot[1] - positions[index + 1][1]
            dy = positions[index + 1][0] - pivot[0]
        elif rotation == 'counterclockwise_xy':
            # Counterclockwise rotation in the XY plane (horizontal plane)
            # Adjust the coordinates accordingly for a 3D rotation
            dx = positions[index + 1][1] - pivot[1]
            dy = pivot[0] - positions[index + 1][0]
        elif rotation == 'clockwise_xz':
            # Clockwise rotation in the XZ plane (vertical plane, front view)
            # Adjust the coordinates accordingly for a 3D rotation
            dx = pivot[2] - positions[index + 1][2]
            dz = positions[index + 1][0] - pivot[0]
        elif rotation == 'counterclockwise_xz':
            # Counterclockwise rotation in the XZ plane (vertical plane, front view)
            # Adjust the coordinates accordingly for a 3D rotation
            dx = positions[index + 1][2] - pivot[2]
            dz = pivot[0] - positions[index + 1][0]
        elif rotation == 'clockwise_yz':
            # Clockwise rotation in the YZ plane (vertical plane, side view)
            # Adjust the coordinates accordingly for a 3D rotation
            dy = pivot[2] - positions[index + 1][2]
            dz = positions[index + 1][1] - pivot[1]
        elif rotation == 'counterclockwise_yz':
            # Counterclockwise rotation in the YZ plane (vertical plane, side view)
            # Adjust the coordinates accordingly for a 3D rotation
            dy = positions[index + 1][2] - pivot[2]
            dz = pivot[1] - positions[index + 1][1]
        else:
            return False  # Invalid rotation

        new_positions[index + 1] = (positions[index + 1][0] + dx, positions[index + 1][1] + dy, positions[index + 1][2] + dz)

        return cstate.is_valid_configuration(new_positions)
    
class State:
    def __init__(self, state):
        self.positions = [(x, y,z) for x, y, z,_ in state]
        self.sequence = [d for _, _,_, d in state]
        self.state = state
        self.energy_matrix = {
                                'HH': -1, 'CC': -5, 'CH': -1, 'HC':-1, 
                                'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
                            }
    def apply_action(self, index, direction):
        current_state = [list(sequence) for sequence in self.state]  # Convert tuples to lists for mutability
        positions = self.positions

        # Check if the action is valid; if not, return the original state immediately
        if index <= 0 or index >= len(positions) - 1:
            return current_state

        pivot = positions[index]
        new_positions = positions.copy()

        dx, dy, dz = 0, 0, 0

        if direction == 'up':
            dy = 1
        elif direction == 'down':
            dy = -1
        elif direction == 'left':
            dx = -1
        elif direction == 'right':
            dx = 1
        elif direction == 'front':
            dz = 1
        elif direction == 'behind':
            dz = -1

        # Calculate the new position based on the previous one
        new_x = pivot[0] + dx
        new_y = pivot[1] + dy
        new_z = pivot[2] + dz

        # Check if the new position is within one unit distance of the previous one
        if abs(new_x - positions[index + 1][0]) <= 1 and abs(new_y - positions[index + 1][1]) <= 1 and abs(new_z - positions[index + 1][2]) <= 1:
            new_positions[index + 1] = (new_x, new_y, new_z)

            # Check if the new configuration is valid
            if not self.is_valid_configuration(new_positions):
                return current_state

            # Update the state with new positions
            for i, position in enumerate(new_positions):
                current_state[i][0] = position[0]
                current_state[i][1] = position[1]
                current_state[i][2] = position[2]

        return current_state

    
    def get_positions(self):
        return [(x, y, z) for x, y, z, _ in self.state]
    
    @staticmethod
    def is_valid_configuration(positions): 
        return len(positions) == len(set(positions))
    
    def reward_function(self):
        energy = 0
        for i in range(len(self.sequence)):
            for j in range(i + 2, len(self.sequence)):  
                distance = abs(self.positions[i][0] - self.positions[j][0]) + abs(self.positions[i][1] - self.positions[j][1]) + abs(self.positions[i][2] - self.positions[j][2])
                if distance == 1:
                    pair = ''.join(sorted([self.sequence[i], self.sequence[j]]))
                    energy += self.energy_matrix[pair]
        return energy

        

    def compute_reward(self, new_state):
        # Save the current state
        current_positions, current_sequence = self.positions, self.sequence

        # Update state to the new state
        self.positions = [(x, y, z) for x, y, z,_ in new_state]
        self.sequence = [z for _, _, _, z in new_state]

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
       # altscore = altheuristic(sequence, action[0])
       # ptscore = currentletter(sequence, action[0])
       # compactscore = compactness_heuristic(state, nextstate)
       # patternscore = folding_heuristic(sequence, action[0])
        rewardscore = cstate.compute_reward(nextstate)
       #  distance = distance_heuristic(state, nextstate)

        #heuristic_scores = [altscore, ptscore, compactscore, patternscore, rewardscore, distance]
        heuristic_scores = rewardscore
        # Calculate the total score using policy weights
        #total = sum(weight * score for weight, score in zip(policy_weights, heuristic_scores))
        new_state = State(nextstate)
        total = -2*new_state.reward_function()
        # Create a vector indicating the sign of each heuristic score
        #sign_vector = [1 if score > 0 else -1 if score < 0 else 0 for score in heuristic_scores]
        sign_vector = 1
        return total, sign_vector



#update policy weight
sequence = "HPHPPHHPHPPHPHHPPHPH"
# policy_weights = [0.1, 0.1, 0.1, 0.1, 0, 0.1]
policy_weights = [ 0.1]
simulator = ProteinFoldingSimulator(sequence, policy_weights)
simulator.run(num_episodes=55, learning_rate=0.01)

print(simulator.bestscore)
print(simulator.beststate)
simulator.visualize_protein()
simulator.plot_weights_and_bond_scores()

