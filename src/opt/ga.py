import os
import sys
from typing import List, Tuple, Set, Dict

current_dir = os.path.dirname(os.path.realpath(__file__))
charac_dir = os.path.join(current_dir, '..', 'charac')
penalty_dir = os.path.join(current_dir, '..', 'penalties')

sys.path.append(current_dir)
sys.path.append(penalty_dir)
sys.path.append(charac_dir)

import mgraph
import boxing
import penalty
import initpop

class individual:
    def __init__(self, initial_solution: List[int], no_genes: int, graph: mgraph.mgraph,
                 subgraph: mgraph.subgraph, target_fsize: int, boxarray: boxing.mbox_array) -> None:
        
        self.m_target_size: int = target_fsize
        self.m_no_genes: int = no_genes

        # the following are references to pre-existing objects, these are not modified
        self.m_graph: mgraph.mgraph = graph
        self.m_subgraph: mgraph.subgraph = subgraph
        self.m_boxarray: boxing.mbox_array = boxarray

        # copy solution
        self.solution: List[int] = initial_solution.copy()
    
    def calculate_score(self):
        return penalty.calculate_score(self)


class population:
    def __init__(self, graph: mgraph.mgraph, subgraph: mgraph.subgraph, target_fsize: int,
                 boxarray: boxing.mbox_array) -> None:
        
        self.m_target_size: int = target_fsize
        self.m_generation: int = 0
        self.m_min_fitness_iters: int = 0
        self.m_best_sol: List[int] = []

        # references to pre-existing objects
        self.m_graph: mgraph.mgraph = graph
        self.m_subgraph: mgraph.subgraph = subgraph
        self.m_boxarray: boxing.mbox_array = boxarray

        no_genes: int = len(subgraph.m_feasible_edges)
        self.m_off_lim_edges: List[int] = [0 for _ in range(no_genes)]
        self.m_best_solution: List[int] = [0 for _ in range(no_genes)]

        initial_population: List[List[int]] = initpop.get_initial_population(graph, subgraph, target_fsize)
        print(f"initial_population size: {len(initial_population)}")
        # self.num_parents_mating: int




        