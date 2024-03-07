import os
import sys
import math
from typing import List, Tuple, Set, Dict

current_dir = os.path.dirname(os.path.realpath(__file__))
charac_dir = os.path.join(current_dir, '..', 'charac')

sys.path.append(current_dir)
sys.path.append(charac_dir)

import mgraph
import param

def eg_vdw(iatom: int, jatom: int, graph: mgraph.mgraph, istatus: bool = False, jstatus: bool = False) -> float:
    atom_types: List[int] = graph.m_atom_types
    
    iatom_type: int = 0 if istatus else atom_types[iatom]
    jatom_type: int = 0 if jstatus else atom_types[jatom]
    Ri: float = param.x1_array[iatom_type]
    Rj: float = param.x1_array[jatom_type]

    ki: float = param.D1sqrt_array[iatom_type]
    kj: float = param.D1sqrt_array[jatom_type]

    kij: float = param.KCAL_TO_KJ * ki * kj
    kaSquared: float = Ri * Rj
    rijSquared: float = graph.distance2(iatom, jatom)

    term6: float = kaSquared / rijSquared
    term6squared: float = term6 * term6 
    term6_final: float = term6squared * term6

    term12: float = term6_final * term6_final
    energy = kij * ((term12) - (2.0 * term6_final))
    return energy


def calculate_evdw_sg(graph: mgraph.mgraph, subgraph: mgraph.subgraph = None) -> None:
    if not subgraph:
        subgraph: mgraph.subgraph = graph.m_subgraphs[0]
    natoms: int = subgraph.m_natoms
    natom_pairs: int = int((natoms * (natoms - 1))/ 2)

    evdw_total: float = 0.0
    for k in range(0, natom_pairs):
        jatom: int = int(natoms - 2 - math.floor(math.sqrt(-8*k + 4*natoms*(natoms-1)-7)/2.0 - 0.5))
        iatom: int = int(k + jatom + 1 - natoms*(natoms-1)/2 + (natoms-jatom)*((natoms-jatom)-1)/2)

        global_jatom: int = subgraph.get_node(jatom)
        global_iatom: int = subgraph.get_node(iatom)

        if graph.distance2(global_iatom, global_jatom) > 144.0:
            continue
        if not graph.is_one_two_bonds_apart(global_iatom, global_jatom):
            

            evdw_total += eg_vdw(global_iatom, global_jatom, graph)

    subgraph.m_energy = evdw_total