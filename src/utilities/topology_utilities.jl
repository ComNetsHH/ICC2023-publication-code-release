#
#  Copyright (c) 2023 Institute of Communication Networks (ComNets),
#                     Hamburg University of Technology (TUHH),
#                     https://www.tuhh.de/comnets
#  Copyright (c) 2023 Leonard Fisser <leonard.fisser@tuhh.de>
# 
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
# 
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
# 
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>.
# 

using Graphs
using Combinatorics

function is_connected_dominating_set(full_graph, mcds_nodes, mcds_candidate_graph)
    N = nv(full_graph)
    nodes_covered = size(unique(vcat(mcds_nodes,reduce(vcat, [neighbors(full_graph,mcds_node) for mcds_node in mcds_nodes]))),1)
    if nodes_covered == N && is_connected(mcds_candidate_graph)
        return true
    else 
        return false
    end
end

function node_set_to_subgraph(connectivity, node_subset)
    node_subset = collect(node_subset)
    d = Array(zeros(size(node_subset,1),size(node_subset,1)))
    for (i_id,i) in enumerate(node_subset)
        for (j_id,j) in enumerate(node_subset)
            d[j_id,i_id] = connectivity[j,i]
            d[i_id,j_id] = connectivity[j,i]
        end
    end
    return Graph(d)
end

function get_node_types(connectivity)
    N = size(connectivity, 1)
    if N <= 16
        return get_node_types_mcds(connectivity) # determine leaves, pseudo-leaves (Farazi paper) and non-leaves, as well as MCDS size
    else 
        return get_node_types_reduced(connectivity) # only finds the "true leaves" and consideres all other nodes non-leaves, does not calculate mcds
    end
end

function get_node_types_reduced(connectivity)
    non_leaves = []
    true_leaves = []
    pseudo_leaves = []
    mcds = []
    best_cds_size = 0
    graph = Graph(connectivity)
    N = nv(graph)
    nodes = 1:N
    for n in nodes
        if sum(connectivity[n,:]) == 1
            append!(true_leaves, n)
        else
            append!(non_leaves, n)
        end
    end
    return best_cds_size, non_leaves, true_leaves, pseudo_leaves, mcds
end
    
function get_node_types_mcds(connectivity)
    non_leaves = []
    true_leaves = []
    pseudo_leaves = []

    graph = Graph(connectivity)
    N = nv(graph)
    nodes = 1:N

    mcds_candidates = combinations(nodes)
    best_cds_size = N
    mcds = []
    # Find a single minimum connected dominating set
    for mcds_candidate in mcds_candidates
        mcds_candidate = collect(mcds_candidate)
        mcds_candidate_graph = node_set_to_subgraph(connectivity, mcds_candidate)
        if is_connected_dominating_set(graph, mcds_candidate, mcds_candidate_graph)
            if length(mcds_candidate) < best_cds_size
                best_cds_size = length(mcds_candidate) 
            end
        end
    end

    # Check all CDS, check if MCDS, union mcds nodes
    for mcds_candidate in mcds_candidates
        mcds_candidate = collect(mcds_candidate)
        mcds_candidate_graph = node_set_to_subgraph(connectivity, mcds_candidate)
        if is_connected_dominating_set(graph, mcds_candidate, mcds_candidate_graph)
            if length(mcds_candidate) == best_cds_size
                union!(non_leaves, collect(mcds_candidate))
                mcds = deepcopy(mcds_candidate)
            end
        end
    end
    true_leaves = findall(vec(sum(connectivity, dims=1).==1))
    pseudo_leaves = setdiff(1:N, true_leaves, non_leaves)
    return best_cds_size, non_leaves, true_leaves, pseudo_leaves, mcds
end
