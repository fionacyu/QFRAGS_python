import os
import sys
from typing import List, Tuple, Set, Dict

current_dir = os.path.dirname(os.path.realpath(__file__))
charac_dir = os.path.join(current_dir, '..', 'charac')
gtraverse_dir = os.path.join(current_dir, '..', 'gtraverse')
sys.path.append(charac_dir)
sys.path.append(current_dir)
sys.path.append(gtraverse_dir)

import mgraph
import gtraverse
import molecule

def calculate_score(self) -> float:
    # member function of inidividual
    natoms: int = subgraph_copy.m_natoms
    global_natoms: int = self.m_graph.m_natoms
    subgraph_copy: mgraph.subgraph = mgraph.subgraph(self.m_subgraph)

    edges: List[Tuple[int, int]] = self.m_graph.m_edges
    feasible_edges: List[int] = self.m_subgraph.m_feasible_edges
    node_sg_nidx: List[int] = self.m_graph.m_node_sg_nidx
    broken_bonds: List[int] = []

    conjugated_systems: Set[int] = set()
    hyperconjugated_systems: Set[int] = set()

    for isol, sol in enumerate(self.solution):
        if sol == 1:
            # edge has been cut
            edge: Tuple[int, int] = edges[feasible_edges[isol]]
            node1: int = edge[0]
            node2: int = edge[1]

            node1_idx: int = node_sg_nidx[node1]
            node2_idx: int = node_sg_nidx[node2]

            subgraph_copy.delete_edge(node1_idx, node2_idx)
            broken_bonds.append(node1_idx)
            broken_bonds.append(node2_idx)

            box_id1: int = self.m_subgraph.get_boxid(node1)
            box_id2: int = self.m_subgraph.get_boxid(node2)

            conjugated_systems1: Set[int] = self.m_boxarray.get_boxes()[box_id1].m_conj_systems
            conjugated_systems2: Set[int] = self.m_boxarray.get_boxes()[box_id2].m_conj_systems

            hyperconjugated_systems1: Set[int] = self.m_boxarray.get_boxes()[box_id1].m_hyper_systems
            hyperconjugated_systems2: Set[int] = self.m_boxarray.get_boxes()[box_id2].m_hyper_systems

            for conj_system in self.m_subgraph.m_conjugated_systems:
                if conj_system in conjugated_systems1 or conj_system in conjugated_systems2:
                    conjugated_systems.add(conj_system)
            
            for hyperconj_system in self.m_subgraph.m_hyperconjugated_systems:
                if hyperconj_system in hyperconjugated_systems1 or hyperconj_system in hyperconjugated_systems2:
                    hyperconjugated_systems.add(hyperconj_system)
    
    fragid: List[int] = [0 for _ in natoms]
    nfrags: int = gtraverse.determine_fragid(fragid, subgraph_copy)

    # define monomers and dimers
    monomer_list: List[molecule.Molecule] = []
    dimer_list: List[molecule.Molecule] = []

    for ifrag in range(0, nfrags):
        monomer: molecule.Molecule = molecule.Molecule([ifrag])
        monomer_list.append(monomer)
    
    for iatom in range(0, natoms):
        fid: int = fragid[iatom]
        atom_idx: int = self.m_subgraph.m_nodes[iatom]
        monomer_list[fid - 1].m_atom_ids.append(atom_idx)
    
    for ifrag in range(0, nfrags): 
        monomer_list[ifrag].set_natoms_nohcap()
    
    # adding hydrogen caps to monomers
    for isol, sol in enumerate(self.solution):
        if sol == 1:
            edge: Tuple[int, int] = edges[feasible_edges[isol]]
            node1: int = edge[0]
            node2: int = edge[1]

            node1_idx: int = node_sg_nidx[node1]
            node2_idx: int = node_sg_nidx[node2]

            fid1: int = fragid[node1_idx]
            fid2: int = fragid[node2_idx]

            if (fid1 != fid2):
                monomer_list[fid1 - 1].m_atom_ids.append(node2 + global_natoms)
                monomer_list[fid2 - 1].m_atom_ids.append(node1 + global_natoms)
            else:
                self.solution[isol] = 0
    
    for ifrag in range(0, nfrags): 
        monomer_list[ifrag].set_natoms()

    # now create dimers
    for imon in range(0, nfrags):
        for jmon in range(0, imon):
            dimer: molecule.Molecule = molecule.Molecule([jmon, imon], monomer_list[jmon], monomer_list[imon],
                                     fragid, self.m_graph)
            dimer_list.append(dimer)

    dimer_energy: List[float] = [0.0 for _ in range(dimer_list)]

    return 0.0