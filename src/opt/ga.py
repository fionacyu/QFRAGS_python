import os
import sys
import math
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

        self.p_comp: float
        self.p_vol: float
        self.p_vrange: float
        self.p_pe: float
        self.p_conj: float
        self.p_hyper: float
    
    def calculate_score(self):
        return penalty.calculate_score(self)
    
    def vol_penalty(self, subgraph_copy: mgraph.subgraph):
        return penalty.vol_penalty(self, subgraph_copy)
    
    def calculate_fragment_volume(self, subgraph_copy: mgraph.subgraph, seed_node: int, colors: List[int]):
        return penalty.calculate_fragment_volume(self, subgraph_copy, seed_node, colors)



class population:
    def __init__(self, graph: mgraph.mgraph, subgraph: mgraph.subgraph, target_fsize: int,
                 boxarray: boxing.mbox_array) -> None:
        
        self.m_target_size: int = target_fsize
        self.m_generation: int = 0
        self.m_min_fitness_iters: int = 0
        self.m_best_sol: List[int] = []
        self.m_ngenes: int = len(subgraph.m_feasible_edges)

        # references to pre-existing objects
        self.m_graph: mgraph.mgraph = graph
        self.m_subgraph: mgraph.subgraph = subgraph
        self.m_boxarray: boxing.mbox_array = boxarray

        self.m_off_lim_edges: List[int] = [0 for _ in range(self.m_ngenes)]
        self.m_best_solution: List[int] = [0 for _ in range(self.m_ngenes)]

        initial_population: List[List[int]] = initpop.get_initial_population(graph, subgraph, target_fsize)
        print(f"initial_population size: {len(initial_population)}")
        
        init_population_size: int = len(initial_population)
        self.num_parents_mating: int

        if (init_population_size > 1 and init_population_size<= 8):
            self.m_num_parents_mating = 2
        elif (init_population_size == 1):
            self.m_num_parents_mating = 1
        elif (init_population_size > 8):
            self.m_num_parents_mating = math.floor(0.25 * init_population_size)
        
        self.m_num_offspring: int = init_population_size - self.m_num_parents_mating
        self.m_mutation_num_genes: int = math.floor(0.1 * self.m_ngenes)
        if (self.m_mutation_num_genes < 1):
            self.m_mutation_num_genes = 0
        
        self.m_solutions: List[individual] = []
        for solution in initial_population:
            indv: individual = individual(solution, self.m_ngenes, self.m_graph, self.m_subgraph, self.m_target_size, self.m_boxarray)
            indv.calculate_score()



        