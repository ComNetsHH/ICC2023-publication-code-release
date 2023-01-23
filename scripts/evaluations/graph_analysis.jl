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

using Graphs, GraphIO, Distributed, DataFrames, Statistics
include("$(srcdir())/solvers/a_star.jl")
include("$(srcdir())/solvers/incremental_a_star.jl")
include("$(srcdir())/solvers/nc_capacities.jl")
include("$(srcdir())/utilities/topology_utilities.jl")

function run_experiment(graph_size, config)
    # Load experiment config
    graphs = readlines("$(datadir())/exp_res/graph$(graph_size)c$(config["cutoff"]).g6")
    graph = GraphIO.Graph6._g6StringToGraph(graphs[config["graph_id"]])
    connectivity = collect(adjacency_matrix(graph))
    solution_knowledge = Dict([("dissemination_delays"=>zeros(graph_size,graph_size)),("t_star",0)]) # Lower Bound Initialisation

    # Prepare solver input, with or without network coding
    if config["max_coding_degree"] == 1
        mcds_cardinality, non_leaves, true_leaves, pseudo_leaves, mcds = get_node_types(connectivity)
        # t_star taken from Farazi paper
        solution_knowledge["t_star"] = graph_size*mcds_cardinality+(graph_size-length(non_leaves))
    else
        # t_star taken from Widmer paper
        solution_knowledge["t_star"] = sum(min_schedule_length_nc(connectivity))
        # We might need to have information from the non-network coding solution, so try to load it
        sys_config = deepcopy(config)
        sys_config["max_coding_degree"] = 1
        sys_results = try wload(datadir("exp_raw/$(graph_size)/$(savename(sys_config,"jld2"))"))
        catch e
            sys_results = missing
        end
        if !ismissing(sys_results)
            first_encodings, first_decodings = get_coding_events(connectivity, sys_results["schedule"])
            solution_knowledge["dissemination_delays"] = first_decodings .- first_encodings # Tighten lower bounds to improve convergence
        end
    end

    # Perform A-Star or Incremental A-Star (Heuristic)
    if (config["max_coding_degree"] > 1) && !config["heuristic"] 
        processing_time = @elapsed best_solution = a_star_search(connectivity, config["max_coding_degree"], config["max_coding_degree"], solution_knowledge)
    else
        processing_time = @elapsed best_solution = incremental_a_star_search(connectivity, config["max_coding_degree"], config["max_coding_degree"],  deepcopy(solution_knowledge))
    end

    schedule = best_solution.schedule
    aoi_components = get_aoi_components(connectivity, schedule)
    
    # Collecting Result Dataframe
    graph_id            = config["graph_id"]
    avg_outdegree       = mean(outdegree(graph))
    graph               = graphs[config["graph_id"]]
    cutoff              = config["cutoff"]
    t                   = length(schedule)
    t_star              = solution_knowledge["t_star"]
    max_coding_degree   = config["max_coding_degree"]
    max_inflight        = config["max_coding_degree"]
    heuristic           = config["heuristic"]

    result = @strdict schedule t t_star aoi_components max_coding_degree heuristic max_inflight graph_id graph processing_time avg_outdegree cutoff
    wsave(datadir("exp_raw/$(graph_size)",savename(config,"jld2")), result)

    @info "Results \n ## Graph ID: $graph_id \n ## Schedule Length: $(t) \n ## Network Coding: $(max_coding_degree>1) \n ## Avg. AoI: $(aoi_components[1])"
    println("")
    return nothing
end