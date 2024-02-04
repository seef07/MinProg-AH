import random
import math
from heuristics.heuristics import altheuristic, currentletter, compactness_heuristic, folding_heuristic, distance_heuristic, hydrophobic_compactness_heuristic, cytosine_compactness_heuristic
import matplotlib.pyplot as plt
import numpy as np 

EPISODE = 1

class ProteinFoldingSimulator:
    """
    A simulator for protein folding using reinforcement learning to optimize folding strategies based on a policy gradient method.
    
    Attributes:
        sequence (str): The sequence of amino acids in the protein.
        energy_matrix (dict): A dictionary mapping pairs of amino acids to their interaction energy.
        policy_weights (list): A list of weights for different heuristics in the policy.
        action_selector (ActionSelector): An instance of ActionSelector to choose actions based on the policy.
        state (list): The current state of the protein, represented as coordinates and amino acid types.
        scores_per_iteration (list): A history of scores for each iteration.
        beststate (list): The best state found during simulation.
        bestscore (float): The score of the best state.
        weights_history (list): A history of policy weights over iterations.
        bond_scores_history (list): A history of bond scores over iterations.
        bestweight (list): The best weight configuration found.
    """
    def __init__(self, sequence, policy_weights):
        """
        Initializes the ProteinFoldingSimulator with a given sequence and initial policy weights.
        """
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
        self.bestweight = None

    def _initialize_state(self):
        return [(i, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]

    def run(self, num_episodes, learning_rate):
        """
        Runs the simulation for a specified number of episodes with a given learning rate.
        
        Parameters:
            num_episodes (int): The number of episodes to run the simulation for.
            learning_rate (float): The learning rate for updating policy weights.
        """
        print("Start run")
        global iteration
        global EPISODE
        randomize_weights_interval = 20
        for episode in range(num_episodes):
            EPISODE = episode
            iteration = 0
            if episode == 1 or episode == 2 or episode == 0:
                EPISODE =1000
                self._run_episode(learning_rate, 1000)
            else:
                self._run_episode(learning_rate, episode)

            #if randomize_weights_interval and episode % randomize_weights_interval == 0 and episode != 0:
             #   self.randomize_weights()


    def _run_episode(self, learning_rate, episode):
        """
        Runs a single episode of the simulation.
        
        Parameters:
            learning_rate (float): The learning rate for updating policy weights.
            episode (int): The current episode number.
        """
        current_state = [(i, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]
        episode_data = []
        for step in range(episode):
            action, heuristic_influence = self.action_selector.select_action(current_state)
            state_instance = State(current_state)
            new_state = state_instance.apply_action(action[0], action[1])  # Apply action
            reward = state_instance.compute_reward(new_state)  # Compute reward
            bondscore = state_instance.reward_function()
            self.bond_scores_history.append(bondscore)
            self.weights_history.append(self.policy_weights.copy())
            if self.bestscore > bondscore:
                self.bestscore = bondscore
                self.beststate = new_state
                self.bestweight = self.policy_weights
            
            episode_data.append((current_state, action, reward, heuristic_influence))
            current_state = new_state
        print("###############################")
        print(self.bestscore)
        print(self.policy_weights)
        print(self.beststate)
        print("###############################")
        self.policy_weights = self.update_policy_weights(self.policy_weights, episode_data, learning_rate)


    def bigrun(self, episode):
        """
        Performs a run using the best weights found from previous simulations.
        
        Parameters:
            episode (int): The number of steps to run in this big run.
        """
        global EPISODE
        EPISODE = episode
        action_class = ActionSelector(self.bestweight)
        current_state = [(i, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]
        for step in range(episode):
            action, heuristic_influence = action_class.select_action(current_state)
            state_instance = State(current_state)
            new_state = state_instance.apply_action(action[0], action[1])  # Apply action
            reward = state_instance.compute_reward(new_state)  # Compute reward
            bondscore = state_instance.reward_function()
            new_state_instance = State(new_state)
            new_bondscore = new_state_instance.reward_function()
            print(bondscore)
            if new_bondscore < bondscore:
                current_state = new_state
            if random.random() < 0.5:
                current_state = new_state

    def randomize_weights(self):
        """
        Randomizes the policy weights by adding a small random value to each.
        """
        self.policy_weights = [w + np.random.uniform(-1, 1) for w in self.policy_weights]

    @staticmethod
    def update_policy_weights(current_weights, episode_data, learning_rate):
        """
        Updates the policy weights based on the episode data and learning rate.
        
        Parameters:
            current_weights (list): The current set of policy weights.
            episode_data (list): The data collected from the current episode.
            learning_rate (float): The learning rate for updating weights.
        
        Returns:
            list: The updated set of policy weights.
        """
        updated_weights = current_weights.copy()

        for state, action, reward, heuristic_influence in episode_data:
            for i in range(len(current_weights)):
                # Adjust weights based on reward and heuristic influence
                if reward > 0 and heuristic_influence[i] > 0 or reward < 0 and heuristic_influence[i] < 0:
                    updated_weights[i] *= (1+learning_rate/2) 
                elif reward > 0 and heuristic_influence[i] < 0 or reward < 0 and heuristic_influence[i] > 0:
                    updated_weights[i] *= (1-learning_rate/4)
               # elif reward == 0 and heuristic_influence[i] > 0 or reward == 0 and heuristic_influence[i] > 0:
                #    updated_weights[i] *= (1-learning_rate/24)
                # Clip the updated weight to ensure it stays within -1 and 1
                updated_weights[i] = np.clip(updated_weights[i], 0.1, 1)
                updated_weights[4] = np.clip(updated_weights[4], 0.4, 1)
        return updated_weights


    def visualize_protein(self):
        """
        Visualizes the protein structure using matplotlib.
        """
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
        """
        Plots the history of policy weights and bond scores over the simulation iterations.
        """
        weights_array = 10*np.array(self.weights_history)
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
    """
    Selects actions based on the current state and policy weights.
    
    Attributes:
        policy_weights (list): The current set of policy weights for decision making.
    """
    
    def __init__(self, policy_weights):
        """
        Initializes the ActionSelector with the given policy weights.
        """
        self.policy_weights = policy_weights

    def select_action(self, state, start_temp=0.8, end_temp=0.01):
        """
        Selects an action for the given state based on the policy weights.
        
        Parameters:
            state (list): The current state of the protein.
            start_temp (float): The starting temperature for simulated annealing.
            end_temp (float): The ending temperature for simulated annealing.
        
        Returns:
            tuple: The selected action and its heuristic influence.
        """
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
        list = []
        for i, part  in enumerate(state):
            actiontrue = (i, True)
            actionfalse = (i, False)
            if self._validate_action(state, actiontrue):
                list.append(actiontrue)
            if self._validate_action(state, actionfalse):
                list.append(actionfalse)
        return list
    
    def _validate_action(self, state, action, clockwise=True): 
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
    """
    Represents the state of the protein during simulation, including its positions and sequence.
    
    Attributes:
        positions (list): The positions of amino acids in the protein.
        sequence (list): The sequence of amino acids in the protein.
        state (list): The combined state information, including positions and types.
        energy_matrix (dict): A dictionary mapping pairs of amino acids to their interaction energy.
    """
    def __init__(self, state):
        """
        Initializes the State with the given state information.
        """
        self.positions = [(x, y) for x, y, _ in state]
        self.sequence = [z for _, _, z in state]
        self.state = state
        self.energy_matrix = {
                                'HH': -1, 'CC': -5, 'CH': -1, 'HC':-1, 
                                'HP': 0, 'PH': 0, 'PP': 0, 'PC': 0, 'CP': 0
                            }
    def apply_action(self, index, clockwise=True):
        """
        Applies an action (rotation) to the protein state.
        
        Parameters:
            index (int): The index of the amino acid around which to rotate.
            clockwise (bool): Whether to rotate clockwise or counterclockwise.
        
        Returns:
            list: The new state after applying the action.
        """
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
        """
        Computes the reward for transitioning to a new state.
        
        Parameters:
            new_state (list): The state after applying an action.
        
        Returns:
            float: The reward for the transition.
        """
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
    """
    Contains methods for evaluating states and actions using various heuristics.
    """
    def compute_heuristic_score(state, action, policy_weights):
        """
        Computes the score for a given action based on heuristics and policy weights.
        
        Parameters:
            state (list): The current state of the protein.
            action (tuple): The action being evaluated.
            policy_weights (list): The current set of policy weights.
        
        Returns:
            float: The total score for the action.
            list: A vector indicating the sign of each heuristic score.
        """
        cstate = State(state)
        nextstate = cstate.apply_action(action[0], action[1])
        sequence = cstate.sequence

        # Compute individual heuristic scores
        altscore = altheuristic(sequence, action[0])
        ptscore = currentletter(sequence, action[0])
        compactscore = compactness_heuristic(state, nextstate)
        patternscore = folding_heuristic(sequence, action[0])
      #  rewardscore = cstate.compute_reward(nextstate)
        distance =15* distance_heuristic(state, nextstate)
        hcompact = 8*hydrophobic_compactness_heuristic(state, nextstate)
        heuristic_scores = [altscore, ptscore, compactscore, hcompact, distance]
        # Calculate the total score using policy weights
        total = sum(weight * score for weight, score in zip(policy_weights, heuristic_scores))
        new_state = State(nextstate)
        total -= 4*new_state.reward_function()
        # Create a vector indicating the sign of each heuristic score
        sign_vector = [1 if score > 0 else -1 if score < 0 else 0 for score in heuristic_scores]

        return total, sign_vector
    

#update policy weight
sequence = "HPHPPHHPHPPHPHHPPHPH"
policy_weights = [1, 1, 1, 0.5, 1.1]

simulator = ProteinFoldingSimulator(sequence, policy_weights)
simulator.run(num_episodes=500, learning_rate=0.01)

print(simulator.bestscore)
print(simulator.beststate)
simulator.visualize_protein()
simulator.plot_weights_and_bond_scores()

