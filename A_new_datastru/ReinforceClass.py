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
        pass

    def visualize(self):
        pass

class State:
    def __init__(self, positions, sequence):
        self.positions = positions
        self.sequence = sequence
