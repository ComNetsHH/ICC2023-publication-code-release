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

using DrWatson
quickactivate("age-optimal-multisource-flooding")

using Graphs, GraphIO, Distributed
include("$(srcdir())/solvers/a_star.jl")
include("$(srcdir())/solvers/nc_capacities.jl")
include("$(srcdir())/utilities/topology_utilities.jl")
include("$(scriptsdir())/evaluations/graph_analysis.jl")

# Setup
graph_size = 5 # Number of nodes
max_coding_degrees = [1,graph_size] # Systematic Coding: 1, Full Network Coding: graph_size
heuristic = [false] # Enabled: divide-and-conquer method using the max-inflight heuristic, use heuristic for graph_size > 7
generate_random_topologies = false # Exhaustive connected graph sets available fir 3 < until graph_size < 8, use random topologies for larger networks 

# Random Topologies
graph_type = "euclidean_graph" # ["euclidean_graph", "erdos_renyi"]
num_random_graphs = 30 # Number of random graphs to produced
cutoffs = collect(0.2:0.05:1.0) # Graph property: cut-off distance (Eucliden Graphs) or density (Erdos Renyi Graphs)

# Runner Config
force_calculation = true # Overwrite any existing results
tmp_dir = "tmp" # Temporary storage for run configs, delete foldor occasionally to prevent hash collision

## -- ##

# Graph Generator
if generate_random_topologies
    num_graphs = num_random_graphs
    for cutoff in cutoffs
        seed = 1
        graphs = []
        for _ in 1:num_random_graphs
            while true # Generate a new graph until it is connected
                seed += 1 
                if graph_type=="euclidean_graph"
                    graph = euclidean_graph(graph_size, 2; cutoff=cutoff, seed=seed)[1]
                elseif  graph_type=="erdos_renyi"
                    graph = erdos_renyi(graph_size, trunc(Int,cutoff*(graph_size*(graph_size-1))/2), seed=seed)
                else
                    println("Wrong graph type selected, available types are: [euclidean_graph, erdos_renyi]")
                end

                if is_connected(graph)
                    push!(graphs, GraphIO.Graph6._graphToG6String(graph)) 
                    break
                end
            end
        end

        # Save graphs to data set
        open("$(datadir())/exp_res/graph$(graph_size)c$(cutoff).g6","w") do f
            [println(f,g[11:end]) for g in graphs]
        end
    end
else
    num_graphs = size(readlines("$(datadir())/exp_res/graph$(graph_size)c.g6"),1) # Only checks how many graphs 
    global cutoffs
    cutoffs = [""]
end

# Experiment Configs
graph_ids = collect(1:num_graphs)
params = dict_list(Dict(
                    :"graph_id" => graph_ids,
                    :"max_coding_degree" => max_coding_degrees,
                    :"heuristic" => heuristic,
                    :"cutoff" => cutoffs))

# Run (configs are saved so that computation can be manually restarted or passed to a batch system)
res = tmpsave(params, projectdir(tmp_dir))
for r in res
    config = load(projectdir(tmp_dir, r), "params")
    if !isfile(datadir("exp_raw/$(graph_size)", savename(config,"jld2"))) || force_calculation
        @info "Running $config"
        run_experiment(graph_size, config)
    end
end