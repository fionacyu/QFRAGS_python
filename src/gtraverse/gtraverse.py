import os
import sys
import queue
from typing import List, Tuple, Set

current_dir = os.path.dirname(os.path.realpath(__file__))
charac_dir = os.path.join(current_dir, '..', 'charac')
monsorter_dir = os.path.join(current_dir, '..', 'mon_sorter')
sys.path.append(charac_dir)
sys.path.append(current_dir)
sys.path.append(monsorter_dir)

import mgraph
import boxing
import conjugated
import mon_sorter

def bfs_connected_component(graph: mgraph.mgraph, subgraph: mgraph.subgraph, seed_node: int, 
    colors: List[int], comp_label: int) -> None: # breadth first search

    subgraph.m_nodes.append(seed_node)
    graph.m_node_sg[seed_node] = comp_label

    bfs_queue: queue.Queue = queue.Queue()
    bfs_queue.put(seed_node)

    while not bfs_queue.empty():
        node_label: int = bfs_queue.get()
        neighbours: List[int] = graph.get_neighbours(node_label)

        for neighbour in neighbours:
            if colors[neighbour] == 0:
                colors[neighbour] = 1
                bfs_queue.put(neighbour)
                subgraph.m_nodes.append(neighbour)
                graph.m_node_sg[neighbour] = comp_label
                graph.m_node_sg_nidx[neighbour] = subgraph.get_natoms() - 1
        colors[node_label] = 2

    print(f"Subgraph has {subgraph.get_natoms()} atoms.")
    subgraph.set_size()

def identify_connected_components(graph: mgraph.mgraph) -> None:
    natoms: int = graph.m_natoms
    colors: List[int] = [0 for _ in range(natoms)]

    comp_label: int = 0
    for iatom in range(0, natoms):
        if colors[iatom] == 0:
            subgraph: mgraph.subgraph = mgraph.subgraph()
            comp_label += 1
            bfs_connected_component(graph, subgraph, iatom, colors, comp_label)
            subgraph.m_offset = iatom
            graph.m_subgraphs.append(subgraph)
    
    if len(graph.m_subgraphs) > 1:
        raise SystemExit(f"Program only works for one molecule at the moment. You have more than one identified.")

def get_unsaturated_neighbours(graph: mgraph.mgraph, node_label: int) -> Set[int]:
    unsaturated_neighbours: Set[int] = set()
    neighbours: Set[int] = graph.get_neighbours(node_label)

    for neighbour in neighbours:
        coordination_number: int = graph.m_coordination_number[neighbour]
        if coordination_number < 4 and coordination_number > 1:
            unsaturated_neighbours.add(neighbour)
    
    return unsaturated_neighbours

def bfs_conjugated(graph: mgraph.mgraph, seed_node: int, conjsys_count: int, colors: List[int], 
                   boxaray: boxing.mbox_array) -> None:
    conjugated_nodes: Set[int] = set()
    conjugated_edges_hash_values: Set[int] = set()
    bfs_queue: queue.Queue = queue.Queue()

    bfs_queue.put(seed_node)
    conjugated_nodes.add(seed_node)

    while not bfs_queue.empty():
        node_label: int = bfs_queue.get()
        unsat_neighbours = get_unsaturated_neighbours(graph, node_label)

        for neighbour in unsat_neighbours:
            if colors[neighbour] == 0:
                colors[neighbour] = 1
                bfs_queue.put(neighbour)
                conjugated_nodes.add(neighbour)

                edge_hash_value: int
                if (node_label < neighbour):
                    edge_hash_value = node_label + graph.m_natoms * neighbour
                else:
                    edge_hash_value = neighbour + graph.m_natoms * node_label

                conjugated_edges_hash_values.add(edge_hash_value)
        colors[node_label] = 2
    
    box_array: List[boxing.mbox] = boxaray.get_boxes()
    if len(conjugated_nodes) >= 3:
        conjsys_size = len(conjugated_nodes)
        conjugated_system = conjugated.conjsys()

        conjugated_system.set_natoms(conjsys_size)
        conjugated_system.set_nedges(len(conjugated_edges_hash_values))

        for conj_node in conjugated_nodes:
            box_id: int = graph.get_boxid(conj_node)
            box_array[box_id].add_conj_sys(conjsys_count)
            conjugated_system.m_nodes.append(conj_node)
            graph.m_conjugation_status[conj_node] = 1
        
        for conj_edge_hash_value in conjugated_edges_hash_values:
            conjugated_system.m_edges_hash_values.append(conj_edge_hash_value)
        
        graph.m_conjugated_systems.append(conjugated_system)
        conjsys_count += 1

def bfs_tree(graph: mgraph.mgraph, ortho_edges: List[Tuple[int, int]], mst_edges: List[Tuple[int, int]], 
             offset: int, natoms: int) -> None: 
    colors: List[int] = [0 for _ in range(natoms)]

    mst_queue: queue.Queue = queue.Queue() 
    mst_queue.put(offset) # offset also refers to the global index of an atom

    while not mst_queue.empty():
        node_label: int = mst_queue.get()
        neighbours: Set[int] = graph.get_neighbours(node_label)

        for neighbour in neighbours:
            bond: Tuple[int, int]
            if node_label < neighbour:
                bond = (node_label, neighbour)
            else:
                bond = (neighbour, node_label)

            if (colors[neighbour - offset] == 0):
                colors[neighbour - offset] = 1
                mst_queue.put(neighbour)
                mst_edges.append(bond)
            elif colors[neighbour - offset] == 1:
                ortho_edges.append(bond)
        colors[node_label - offset] = 2
            
def shortest_path_bfs(graph: mgraph.mgraph, node1: int, node2: int, offset: int) -> List[int]:
    shortest_path: List[int] = []
    seen: Set[int] = set()
    path_list: List[List[int]] = []

    seen.add(node1 - offset)
    first_path: List[int] = []
    first_path.append(node1 - offset)
    path_list.append(first_path)

    path_idx: int = 0
    while (path_idx < len(path_list)):
        shortest_path = path_list[path_idx]
        last_node: int = shortest_path[-1]
        neighbours: Set[int] = graph.get_neighbours(last_node)

        node2_idx: int = node2 - offset
        if node2_idx in neighbours:
            shortest_path.append(node2_idx)
            return shortest_path
        
        for neighbour in neighbours:
            if neighbour not in seen:
                new_path: List[int] = shortest_path.copy()
                new_path.append(neighbour)
                path_list.append(new_path)
                seen.add(neighbour)
        path_idx += 1

    return shortest_path

def min_cycle(graph: mgraph.mgraph, shortest_path: List[int], nnodes: int, cycle_idx: int,
              boxarray: boxing.mbox_array, subgraph_idx: int) -> None:
    offset: int = graph.m_subgraphs[subgraph_idx].m_offset
    box_array: List[boxing.mbox] = boxarray.get_boxes()

    shortest_path = [x + offset for x in shortest_path]
    
    cutoff: int = nnodes + offset

    for inode in range(0, len(shortest_path)):
        if shortest_path[inode] < cutoff:
            shortest_path[inode] = shortest_path[inode]
        else:
            shortest_path[inode] -= nnodes
    graph.m_subgraphs[subgraph_idx].setup_cycle(cycle_idx, len(shortest_path) - 1)

    boxes: Set[int] = set()
    for inode in range(0, len(shortest_path) - 1):
        node1: int = shortest_path[inode]
        node2: int = shortest_path[inode + 1]

        graph.m_ring_status[node1] = len(shortest_path) - 1
        graph.m_ring_status[node2] = len(shortest_path) - 1

        if (node1 < node2):
            graph.m_subgraphs[subgraph_idx].m_cycles[cycle_idx].m_edges[2 * inode] = node1
            graph.m_subgraphs[subgraph_idx].m_cycles[cycle_idx].m_edges[2 * inode + 1] = node2
        else:
            graph.m_subgraphs[subgraph_idx].m_cycles[cycle_idx].m_edges[2 * inode] = node2
            graph.m_subgraphs[subgraph_idx].m_cycles[cycle_idx].m_edges[2 * inode + 1] = node1

        box_idx: int = graph.get_boxid(node1)
        boxes.add(box_idx)

    graph.m_subgraphs[subgraph_idx].setup_cycle_boxes(cycle_idx, len(boxes))
    for box_id in boxes:
        graph.m_subgraphs[subgraph_idx].m_cycles[cycle_idx].m_boxes.append(box_id)
        box_array[box_id].m_cycles.add(cycle_idx)

def shortest_path_hyperconjugated(node1: int, node2: int, graph: mgraph.mgraph, shortest_path: List[int]) -> None:
    path_dist: int = 0

    if node1 in graph.get_neighbours(node2):
        path_dist = 1
    elif node1 in graph.get_neighbours2(node2):
        path_dist = 2
    elif node1 in graph.get_neighbours3(node2):
        path_dist = 3
    
    node1_neighbours: Set[int] = graph.get_neighbours(node1)
    node2_neighbours: Set[int] = graph.get_neighbours(node2)

    match path_dist:
        case 1:
            shortest_path.append(node1)
            shortest_path.append(node2)
            return

        case 2:
            for neighbour in node1_neighbours:
                if neighbour in node2_neighbours:
                    shortest_path.append(node1)
                    shortest_path.append(neighbour)
                    shortest_path.append(node2)
                    return 
        
        case 3:
            for neighbour1 in node1_neighbours:
                for neighbour2 in node2_neighbours:
                    neighbour2_neighbours: Set[int] = graph.get_neighbours(neighbour2)
                    if neighbour1 in neighbour2_neighbours:
                        shortest_path.append(node1)
                        shortest_path.append(neighbour1)
                        shortest_path.append(neighbour2)
                        shortest_path.append(node2)
                        return

def bfs_sg(subgraph: mgraph.subgraph, seed_node: int, colors: List[int], target_size: int) -> int:
    subgraph_nodes: int = 1

    bfs_queue: queue.Queue = queue.Queue()
    bfs_queue.put(seed_node)

    while not bfs_queue.empty():
        node_label: int = bfs_queue.get()
        neighbours: Set[int] = subgraph.get_neighbours(node_label)
        for neighbour in neighbours:
            if colors[neighbour] == 0:
                colors[neighbour] = 1
                bfs_queue.put(neighbour)
                subgraph_nodes += 1
                if subgraph_nodes / target_size >= 0.6:
                    return subgraph_nodes
        colors[node_label] = 2
    return subgraph_nodes

def bfs_fragid(colors: List[int], fragid: List[int], subgraph_copy: mgraph.subgraph, 
               fragid_idx: int, seed_node: int) -> None:
    bfs_queue: queue.Queue = queue.Queue()
    bfs_queue.put(seed_node)
    fragid[seed_node] = fragid_idx

    while not bfs_queue.empty():
        node_label: int = bfs_queue.get()
        neighbours: Set[int] = subgraph_copy.get_neighbours(node_label)
        for neighbour in neighbours:
            if colors[neighbour] == 0:
                colors[neighbour] = 1
                fragid[neighbour] = fragid_idx
                bfs_queue.put(neighbour)
        colors[node_label] = 2


def determine_fragid(fragid: List[int], subgraph_copy: mgraph.subgraph) -> int:
    natoms: int = subgraph_copy.m_natoms
    colors: List[int] = [0 for _ in range(natoms)]

    fragid_idx: int = 0
    for iatom in range(0, natoms):
        if colors[iatom] == 0:
            fragid_idx += 1
            bfs_fragid(colors, fragid, subgraph_copy, fragid_idx, iatom)

    return fragid_idx # equiv to natoms

def patch_initial_pop_bfs(mon_adjacency: List[Set[int]], seed_mon: int, colors: List[int], frag_idx:int,
                          new_mon_fragids: List[int], new_mon_size: List[int], target_size: int) -> None:
    fsize: int = new_mon_size[seed_mon]
    new_mon_fragids[seed_mon] = frag_idx

    if fsize / target_size >= 0.5:
        return
    
    bfs_queue: queue.Queue = queue.Queue()
    bfs_queue.put(seed_mon)
    while not bfs_queue.empty():
        mon_label: int = bfs_queue.get()
        mon_neighbours: Set[int] = mon_adjacency[mon_label]
        monomer_list: List[mon_sorter.monomer] = []
        for mon_neighbour in mon_neighbours:
            mon: mon_sorter.monomer = mon_sorter.monomer(new_mon_size[mon_neighbour], mon_neighbour)
            # distance is the no. of atoms

            monomer_list.append(mon)
        monomer_list.sort(key=lambda x: x.distance)
        for imon in range(0, len(monomer_list)):
            imon_fsize: int = monomer_list[imon].distance
            if colors[monomer_list[imon].idx] == 0 and (imon_fsize + fsize)/target_size <= 1.5:
                new_mon_fragids[monomer_list[imon].idx] = frag_idx
                fsize += monomer_list[imon].distance
                colors[monomer_list[imon].idx] = 1
                bfs_queue.put(monomer_list[imon].idx)
                
                if fsize/target_size > 0.9:
                    return
        colors[mon_label] = 2
