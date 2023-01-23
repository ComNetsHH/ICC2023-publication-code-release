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

using Graphs, Combinatorics
using JuMP, GLPK, LinearAlgebra
#using Gurobi # Consider using Gurobi for large problems
include("$(srcdir())/utilities/topology_utilities.jl")

function min_schedule_length_nc(connectivity, inactive_nodes = [])
    """Returns the required node capacities for full network coded flooding. Taken from Jörg Widmer paper."""
    model = Model(GLPK.Optimizer)
    g = Graph(connectivity)
    
    N = nv(g)
    @variable(model, C[1:N], Int)
    for c in C
        @constraint(model, c >= 0)
    end

    for cut in combinations(collect(1:N))
        S = cut
        S_prime = setdiff(collect(1:N),S)
        if size(S,1) < N 
            cut_mask = zeros(Int,N)
            cut_mask[S] .= 1
            cut_mask[S_prime] .= 2
            Ns = unique(cat([[n.src n.dst] for n in karger_cut_edges(g, cut_mask)]...,dims=2))
            intersect!(Ns,S)
            setdiff!(S,inactive_nodes)
            @constraint(model,sum(C[Ns]) >= size(S,1))
        end
    end
    @objective(model, Min, sum(C))
    optimize!(model)
    return JuMP.value.(C)
end

function get_dynamic_opt_problem(connectivity)
    """Returns a dynamic optimization problem, for which the set of active nodes can be adjusted."""
    model = Model(GLPK.Optimizer)
    #model = Model(Gurobi.Optimizer) For large problems, Gurobi's performance is required

    set_silent(model)
    g = Graph(connectivity)
    
    N = nv(g)
    model[:C] = @variable(model, C[1:N], Int)
    model[:node_active] = @variable(model, node_active[1:N], Int)

    for c in C
        @constraint(model, c >= 0)
    end

    for cut in combinations(collect(1:N))
        S = cut
        S_prime = setdiff(collect(1:N),S)
        if size(S,1) < N 
            cut_mask = zeros(Int,N)
            cut_mask[S] .= 1
            cut_mask[S_prime] .= 2
            Ns = unique(cat([[n.src n.dst] for n in karger_cut_edges(g, cut_mask)]...,dims=2))
            intersect!(Ns,S)
            @constraint(model,sum(C[Ns]) >= sum(node_active[S]))
        end
    end

    @objective(model, Min, sum(C))
    return model
end