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

class ActionSelector:
    def __init__(self, policy_weights):

        self.policy_weights = policy_weights

    def select_action(self, state):
        pass


class State:
    def __init__(self, state):
        self.positions = [(x, y) for x, y, _ in state]
        self.sequence = [z for _, _, z in state]
        self.state = state

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
    
    def is_valid_configuration(positions): 
        return len(positions) == len(set(positions))

class Heuristic:
    def evaluate(self, state):
        pass


##rewardfunction
    
##compute reward
    
##policy_gradient_main_loop

#extract_positions
    
#select_action
    
#apply action

#validate action

#get possible actions

# is valid

#update policy weight

#compute_heuristic_score

#compute h infl

