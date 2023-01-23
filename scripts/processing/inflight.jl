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

using DataFrames, Statistics, PyPlot, Graphs, GraphIO

include("$(srcdir())/utilities/aoi_utilities.jl")

function pl_inflight(graph_sizes)
    fig, ax = PyPlot.subplots()
    max_inflight_data = []

    labels = collect(1:6)
    width = 0.8  # the width of the bars

    x = collect(1:length(labels))

    for (i,graph_size) in enumerate(graph_sizes)
        df = collect_results(datadir("exp_raw/$graph_size"))
        sys_nc_gdf = groupby(df, :max_inflight, sort=true)
        nc_df = sys_nc_gdf[2]

        g_max_inflight_data = [get_max_inflight_data(get_coding_events(adjacency_matrix(GraphIO.Graph6._g6StringToGraph(nc_df[i,:].graph)),nc_df[i,:].schedule)...) for i in 1:nrow(nc_df)]
        
        push!(max_inflight_data, g_max_inflight_data)
        hist_counts = [length(findall(g_max_inflight_data .== label)) for label in labels]./length(g_max_inflight_data)
        hist_counts[hist_counts .== 0] .= -0.01
        spacing = width/length(graph_sizes)
        ax.bar(x .- spacing*trunc(Int,length(graph_sizes)/2) .+ spacing*i, hist_counts, width=width/length(graph_sizes))

        for x in labels
            print("$(trunc(hist_counts[x]*100, digits=2)) ")
        end
        println(sum(labels[hist_counts.>0].*hist_counts[hist_counts.>0]))
    end
    ax.set_xticks(x .+ (width/length(graph_sizes)/2), labels)
    ax.set_axisbelow(true)
    ax.grid(which="major",  color="k", linestyle="-", linewidth=0.25, axis="y")

    # Set position of bar on X axis
    xlabel("Maximum Inflight")
    ylabel("Relative Frequency")
    println(plotsdir("inflight.png"))
    savefig(plotsdir("inflight.png"))
end

graph_sizes = [4,5,6,7]
pl_inflight(graph_sizes)