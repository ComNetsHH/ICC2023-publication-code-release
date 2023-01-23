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

@everywhere begin
    using DrWatson
    quickactivate("age-optimal-multisource-flooding")

    using DataStructures
    include("$(srcdir())/utilities/aoi_utilities.jl")

    mutable struct IncrementalAStarState
        schedule::Vector{Tuple{Int64, Int64}} # The sequence of actions
        payload_schedule::Vector{Tuple{Int64,Vector{Int8}}} # Also sequence of actions, but with binary payload vector
        decoders::Array{Int8,3} # The decoder matrices of every node
        decoding_levels::Matrix{Int16} # Logs the first time data was decoded at every node -> matrix
        first_encodings::Vector{Int16} # Logs the very first transmissions of data -> vector
        content_purged::BitVector # Data fully disseminated is removed, so we keep track here
        
        f_score::Float32 # a_star variable
        g_score::Float32 # a_star variable
        h::Float32 # a_star variable
        dist::Float32 # a_star variable
        
        end_state::Bool # Schedule end
        present_content::Vector{Int8}
        active_content_map::Dict{Int64, Int64}
    end

    Base.show(io::IO, x::IncrementalAStarState) = print(io, "F: $(x.f_score), G: $(x.g_score), Schedule: $(x.schedule)")
    Base.hash(x::IncrementalAStarState) = hash(deepcopy(x.decoders))
    Base.deepcopy(m::IncrementalAStarState) = IncrementalAStarState(copy(m.schedule),
                                                                    copy(m.payload_schedule),
                                                                    copy(m.decoders),
                                                                    copy(m.decoding_levels),
                                                                    copy(m.first_encodings),
                                                                    copy(m.content_purged),
                                                                    copy(m.f_score),
                                                                    copy(m.g_score),
                                                                    copy(m.h),
                                                                    copy(m.dist),
                                                                    copy(m.end_state),
                                                                    copy(m.present_content),
                                                                    copy(m.active_content_map))
    

    function evaluate_incremental_state!(connectivity, state, solution_knowledge)
        N = size(connectivity,1)
        step = size(state.schedule,1)
        decoders = state.decoders
        good_action = false

        if step != 0
            action = state.schedule[step]
            tx_id = action[1]
            good_action = apply_action!(connectivity, decoders, action)
            if good_action
                reduced_payload = action_to_coded_payload(decoders, action)
                payload = zeros(N)
                for x in state.active_content_map
                    if x[2] != 0
                        if reduced_payload[x[1]] == 1; payload[x[2]] = 1 end
                    end
                end 
                push!(state.payload_schedule,(tx_id, payload))
                if (state.first_encodings[tx_id] == typemax(typeof(state.first_encodings[tx_id]))) && (payload[tx_id] == 1) 
                    state.first_encodings[tx_id] = step
                end
            else
                state.dist = Inf
                return nothing
            end
        end
        
        dist = 0.0 # Cost associated with this action
        h_dist = 0 # Remaining expected costs

        is_decoded = get_is_decoded(decoders, state.active_content_map)
        decoded_indices = findall(is_decoded.==1)
        for index in decoded_indices
            if state.decoding_levels[index[1], index[2]] == 32767 #TODO: improve magic number
                state.decoding_levels[index[1],index[2]] = step + 1
                dist += 1
            end
        end

        before_purge = deepcopy(state.content_purged)
        purge_disseminated_contents!(decoders, state.content_purged, state.active_content_map)
        if (sum(before_purge) < sum(state.content_purged))
            a_p = findall(state.content_purged .!= before_purge)
            for a in a_p
                for x in state.active_content_map
                    if x[2] == a
                        state.active_content_map[x[1]] = 0 
                    end
                end
            end
        end

        dissemination_delays = solution_knowledge["dissemination_delays"]

        for c_id = findall(.!state.content_purged)
            if state.first_encodings[c_id] < 32767 #TODO:
                decoding_status = is_decoded[c_id,:]
                dist += (N-sum(decoding_status)) # Inflight, ages additionaly, for every non-delivered-to node
                # We will at least collect the dissemination_delay, so we should consider this in the heuristic
                for d_id = 1:N
                    if decoding_status[d_id] == 0
                        h_dist += max(0, dissemination_delays[c_id, d_id] - (step-state.first_encodings[c_id]))
                    end
                end
            else
                # we did not yet encode, and we will at least collect the content_ids dissemination delay as dist
                h_dist += sum(dissemination_delays[c_id,:]) - dissemination_delays[c_id,c_id]
            end
        end

        if sum(state.decoders) == 0
            state.end_state = true
            state.h = 0
        else
            min_tx_left = max(1, solution_knowledge["t_star"]-step)
            state.h = (N*(N-1)/2)*min_tx_left + h_dist
        end

        # Add additional aging due to increase in update interval
        dist += (N*(N-1)/2)
        state.dist = dist

        return nothing
    end

    function evaluate_merge(loop_state, merge, action_set, connectivity, max_coding_degree, max_inflight, solution_knowledge, t_star, star_epsilon, true_leaves)
        N = size(connectivity,1)
        reduced_loop_state = deepcopy(loop_state)
        reduced_loop_state.decoders = zeros(Int8, max_inflight, max_inflight, N)
        reduced_loop_state.end_state = false
        solution_knowledge["t_star"] = t_star
        new_content_map = Dict([x=>0 for x in 1:max_inflight])
        union!(reduced_loop_state.present_content,merge)

        # To improve performance and memory footprint, we always just consider a max_inflight x max_inflight decoder matrix
        # This means, that for max_inflight = 2, all decoders are just a 2x2 matrix, which is very performant
        # However, this also means we have to keep track of which entries corrospond to which node's data, eg. index 1 of the matrix represents data from node 7, index 2 from node 1
        # The following code, makes sure, that the mapping is taken care off, although the code is very confuse

        # Transform status of partial disseminated info to reduced form
        for content_id in reduced_loop_state.present_content
            if sum(loop_state.decoders[content_id, :, :]) > 0
                # find unused dict entry
                new_entry = missing
                for x in new_content_map
                    if x[2] == 0
                        new_content_map[x[1]] = content_id
                        new_entry = x[1]
                        break
                    end
                end
                node_ranks = get_node_ranks(reduced_loop_state.decoders)
                for node_id = 1:N
                    row_entry = node_ranks[node_id]+1
                    for row_id = 1:N
                        if loop_state.decoders[content_id, row_id, node_id] == 1
                            reduced_loop_state.decoders[new_entry, row_entry, node_id] = loop_state.decoders[content_id, row_id, node_id]
                            row_entry += 1
                        end
                    end
                end
            end
        end

        reduced_loop_state.active_content_map = deepcopy(new_content_map)

        step = 1
        open_set_pq = PriorityQueue(hash(reduced_loop_state)=>reduced_loop_state.f_score)
        open_set = Dict(hash(reduced_loop_state)=>reduced_loop_state)
        best_solution = missing

        # This is basically the same code as in the a_star version, could be refactored
        while !isempty(open_set_pq)
            if mod(step, 5000) == 0
                println("Heap size: $(length(open_set_pq)), Dict size: $(length(open_set))")
            end

            best = dequeue!(open_set_pq)
            candidate = open_set[best]
            if candidate.end_state
                best_solution = deepcopy(candidate)
                break
            end
            
            branching_actions = valid_actions(action_set, candidate, max_coding_degree, [])
            
            # grow set with all possibilities
            for branching_action in branching_actions
                branch_candidate = deepcopy(candidate)
                branch_candidate.g_score = Inf

                # Add action to schedule
                push!(branch_candidate.schedule, branching_action)
                evaluate_incremental_state!(connectivity, branch_candidate, solution_knowledge)

                if !isinf(branch_candidate.dist)
                    tentative_g_score = candidate.g_score + branch_candidate.dist
                    check = get(open_set, hash(branch_candidate), missing)

                    if ismissing(check)
                        check = branch_candidate
                    end

                    if tentative_g_score < check.g_score
                        check = deepcopy(branch_candidate)
                        check.g_score = tentative_g_score
                        check.f_score = tentative_g_score + star_epsilon*check.h
                        open_set[hash(check)] = check
                        open_set_pq[hash(check)] = check.f_score
                    end
                end
            end
            step += 1
        end
        return best_solution
    end
end

function get_content_costs(connectivity)
    # To evaluate the "usefullness" of merging a specific set of data, a baseline performance for each one of them has to be established.
    # Therefore, a systematic code solution is computed here
    # Searching for an optimal dissemination schedule is very fast, since all decoder matrices are just 1x1's
    N = size(connectivity,1)
    solution_knowledge = Dict([("dissemination_delays"=>zeros(N,N)),("t_star",0)])
    content_costs = zeros(N)
    action_set = vec(collect(Iterators.product(1:N,1)))
    sys_merges = collect(enumerate(collect(1:N)))
    start_decoders = zeros(Int8, N, N, N)
    [start_decoders[n,1,n] = 1 for n = collect(1:N)]
    sys_state = IncrementalAStarState([],[],start_decoders,ones(Int16,N,N).*32767,ones(Int16,N).*32767,falses(N),0.0,0.0,0.0,0.0,false,[],Dict([x=>0 for x in 1:1]))                            
    sys_combination_log = map(merge -> evaluate_merge(deepcopy(sys_state), merge[1], action_set, connectivity, 1, 1, solution_knowledge, 0, 0.0, []), sys_merges)
    for x in sys_combination_log
        content_costs[x.present_content[1]] = x.g_score
    end
    return content_costs
end


function incremental_a_star_search(connectivity, max_coding_degree, max_inflight, solution_knowledge, star_epsilon=0.0)
    N = size(connectivity,1)
    start_decoders = zeros(Int8, N, N, N)
    [start_decoders[n,1,n] = 1 for n = collect(1:N)]
    accumulated_state = IncrementalAStarState([],[],start_decoders,ones(Int16,N,N).*32767,ones(Int16,N).*32767,falses(N),0.0,0.0,0.0,0.0,false,[],Dict([x=>0 for x in 1:max_inflight]))
    content_costs = get_content_costs(connectivity)
    action_set = vec(collect(Iterators.product(1:N,1:(2^max_inflight)-1)))
    true_leaves = []

    # We are keeping track of evaluated merges, to maybe skip computations in the future 
    cached_merges = []
    minimum_purges = 1
    while sum(accumulated_state.content_purged) < N
        # Calculate all possibile content merges
        num_draws = min(max_inflight - (length(accumulated_state.present_content)-sum(accumulated_state.content_purged)), N-length(accumulated_state.present_content))
        draw_from = setdiff(collect(1:N), accumulated_state.present_content)

        merges = []
        if (num_draws == max_inflight) && !isempty(cached_merges)
            for merge in cached_merges
                if !any([m in accumulated_state.present_content for m in merge[2]])
                    merges = [(1,deepcopy(merge[2]))]
                    break
                end
            end
        end
        if isempty(merges)
            merges = collect(enumerate(combinations(draw_from,num_draws)))
        end

        t_stars = zeros(size(merges,1))
        # Calculate the aoi values for each merge
        combination_log = pmap(merge -> evaluate_merge(deepcopy(accumulated_state), merge[2], action_set, connectivity, max_coding_degree, max_inflight, solution_knowledge, t_stars[merge[1]], star_epsilon, true_leaves), merges)
        # Take the merge, which resulted in the best relativ performance compared to the summed systematic solutions of each merge
        sorter = sortperm(combination_log, by = v -> sum(content_costs[v.present_content]) - v.g_score, rev=true)
        combination_log = combination_log[sorter]
        merges = merges[sorter]

        if isempty(cached_merges)
            cached_merges = deepcopy(merges)
        end

        # Rollback back the state to the first instance where one of the merge data was fully disstributed, this is the starting point for all future merges
        best_partial_state = deepcopy(combination_log[1])
        roll_back_state = IncrementalAStarState([],[],zeros(Int8, N, N, N),ones(Int16,N,N).*32767,ones(Int16,N).*32767,falses(N),0.0,0.0,0.0,0.0,false,[],Dict([x=>x for x in 1:N]))
        roll_back_state.present_content = deepcopy(best_partial_state.present_content)
        for x = 1:N
            roll_back_state.decoders[x,1,x] = 1
        end
        evaluate_incremental_state!(connectivity, roll_back_state, solution_knowledge)

        for x in best_partial_state.payload_schedule
            if sum(roll_back_state.content_purged) >= minimum_purges
                break
            end
            roll_back_action = coded_payload_to_action(roll_back_state.decoders, x[1], x[2])
            push!(roll_back_state.schedule, roll_back_action)
            evaluate_incremental_state!(connectivity, roll_back_state, solution_knowledge)
            tentative_g_score = roll_back_state.g_score + roll_back_state.dist
            roll_back_state.g_score = tentative_g_score
            roll_back_state.f_score = tentative_g_score + star_epsilon*roll_back_state.h
        end

        accumulated_state = deepcopy(roll_back_state)
        minimum_purges = sum(roll_back_state.content_purged) + 1
        sleep(1)
    end

    accumulated_state.schedule = payload_schedule_to_action_schedule(connectivity, accumulated_state.payload_schedule)
    return accumulated_state
end