from . import Heuristic

from protein_folding.protein import Protein
from protein_folding.node import Node


# Extra factor to multiply results with so they're not all in range {0.9, 1.]
_MULT_FACTOR = 7


class Potential(Heuristic):
    def __init__(self, *args):
        super().__init__(*args)

    def run(self):
        score = 0
        for i, node in enumerate(self.protein.nodes[:-1]):
            if node.letter == 'P':
                continue
            if node.ghost:
                break

            for other_node in self.protein.nodes[i+1:]:
                if other_node.letter == 'P':
                    continue
                if other_node.ghost:
                    break

                delta_vec = other_node.pos - node.pos
                len_sq = delta_vec.len_sq()
                score += 1 / len_sq

        self.score_per_direction.append(score)

    def interpret(self) -> list[float]:
        """
        Normalises all scores such that max(scores) == 1
        """
        _min = min(self.score_per_direction)
        _max = max(self.score_per_direction)

        if _max == 0:
            return [0. for _ in self.score_per_direction]

        scores_norm = [value / _max for value in self.score_per_direction]  # {, 1]
        scores_inv = [1 - value for value in scores_norm]  # [0, }
        positive_norm_factor = 1 - max(scores_inv)
        scores_corr = [(value + positive_norm_factor) * _MULT_FACTOR - _MULT_FACTOR + 1 for value in scores_inv]  # {, 1]
        return scores_corr

    def _find_sources(self, node_idx: int) -> list[Node]:
        return [node for node in self.protein.nodes[:node_idx] if node.letter != 'P']
    
    
    def calculate_score_for_state(protein_state):
        score = 0
        for i, node in enumerate(protein_state.nodes[:-1]):
            if node.letter == 'P' or node.ghost:
                continue

            for other_node in protein_state.nodes[i+1:]:
                if other_node.letter == 'P' or other_node.ghost:
                    continue

                delta_vec = other_node.pos - node.pos
                len_sq = delta_vec.len_sq()
                if len_sq != 0:
                    score += 1 / len_sq

        return score * _MULT_FACTOR  
