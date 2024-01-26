import random
import math

class ProteinFoldingSimulator:
    def __init__(self, sequence, energy_matrix, policy_weights):
        self.sequence = sequence
        self.energy_matrix = energy_matrix
        self.policy_weights = policy_weights
        self.action_selector = ActionSelector(policy_weights)
        self.state = self._initialize_state()
        self.scores_per_iteration = []

    def _initialize_state(self):
        return [(i, 0, amino_acid) for i, amino_acid in enumerate(self.sequence)]

    def run(self, num_episodes, learning_rate):
        for episode in range(num_episodes):
            self._run_episode(learning_rate)

    def _run_episode(self, learning_rate):
        current_state = self.state
        episode_data = []
        for step in range(10000):
            action = ActionSelector.select_action(self.policy_weights, current_state)  # Use heuristic
            new_state = apply_action(current_state, action[0], action[1])  # Apply action
            reward = compute_reward(new_state, current_state)  # Compute reward
            bondscore = reward_function(extract_positions(new_state), sequence)
            if bestscore > bondscore:
                bestscore = bondscore
                beststate = new_state
            
            episode_data.append((current_state, action, reward))
            current_state = new_state


    def visualize(self):
        pass

iteration = 0
class ActionSelector:
    def __init__(self, policy_weights):
        self.policy_weights = policy_weights

    def select_action(self, state, policy_weights, start_temp=1.0, end_temp=0.01):
        global iteration
        iteration += 1
        temp = start_temp - iteration * (start_temp - end_temp) / 10000             ###UPDATE
        possible_actions = self._get_possible_actions(state)
        action_scores = []

        for action in possible_actions:
            score = compute_heuristic_score(state, action, policy_weights)
            action_scores.append(score)

        max_score = max(action_scores)
        best_actions = [action for action, score in zip(possible_actions,action_scores) if score == max_score]
        selected_action = random.choice(best_actions)
        new_state = state.apply_action(state, selected_action[0], selected_action[1])  
        new_energy = state.reward_function(new_state) 
        current_energy = state.reward_function(state)     
        exp_argument = (new_energy - current_energy) / temp

        capped_argument = min(exp_argument, 1) 

        if random.random() < math.exp(capped_argument):
            selected_action = random.choice(possible_actions)

        return selected_action

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
    
    def _validate_action(state, action, clockwise=True): 
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
    def apply_action(self, index, clockwise = True):
        current_state = [list(sequence) for sequence in self.state]  # Convert tuples to lists for mutability
        positions = self.positions

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

        if self.is_valid_configuration(new_positions):   
            for i, position in enumerate(new_positions):
                current_state[i][0] = position[0]
                current_state[i][1] = position[1]
            return current_state

        return current_state
    
    def get_positions(self):
        return [(x, y) for x, y, _ in self.state]
    
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



##policy_gradient_main_loop

#select_action  ---> Heuristis etcetera



#update policy weight

#compute_heuristic_score

#compute h infl

