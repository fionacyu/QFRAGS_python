import os
import sys
import math
from typing import List, Tuple, Set, Dict

current_dir = os.path.dirname(os.path.realpath(__file__))
charac_dir = os.path.join(current_dir, '..', 'charac')
sys.path.append(charac_dir)
sys.path.append(current_dir)

import mgraph
import boxing
import ga

class Optimiser_Frag:
    def __init__(self, graph: mgraph.mgraph, boxarray: boxing.mbox_array) -> None:
        self.m_target_size: int
        self.m_edge_solution: List[int] = [0 for _ in range(graph.m_nedges)]

        # the following are references to pre-existing objects
        self.m_graph: mgraph.mgraph = graph
        self.m_boxarray: boxing.mbox_array = boxarray
        self.m_starting_subgraph: mgraph.subgraph = graph.m_subgraphs[0]
        self.m_target_size: int = graph.m_target_frag_size
        self.m_subgraphs: List[mgraph.subgraph] = []
    

    def frag(self, subgraph: mgraph.subgraph, child_subgraph: List[mgraph.subgraph], first: bool) -> None:
        natoms: int = subgraph.m_natoms
        size_ratio: float = natoms / self.m_target_size

        target_size: int # we fragment recursively

        if (size_ratio <= 1.0):
            raise SystemExit("Your fragment size is greater than system size.")
        elif (size_ratio > 1.0 and size_ratio <= 1.5):
            target_size = int(natoms / 2.0)
        elif (size_ratio > 1.5 and size_ratio <= 6.0):
            target_size = self.m_target_size
        elif (size_ratio >= 6.0):
            target_size = math.floor(natoms / 5.0)
        
        converged: bool = False
        pop: ga.population = ga.population(self.m_graph, subgraph, target_size, self.m_boxarray)


    def run(self) -> None:
        self.frag(self.m_starting_subgraph, self.m_subgraphs, True)