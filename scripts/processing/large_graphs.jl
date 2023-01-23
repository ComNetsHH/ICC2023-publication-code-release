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

using DataFrames, Statistics, PyPlot, Graphs, GraphIO, Colors

function rounds(value, step)
    return round(round(value / step) * step,digits=2);
end

function pl_graphsize_degree_gain(graph_sizes, group_by, subresult_path)
    x_values = []
    y_values = []
    z_values = []
    sample_values =  []

    ## Set up LaTeX fonts
    plt.rcParams["text.usetex"] = true
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serifx"] =  ["Computer Modern Roman"]
    plt.rcParams["font.size"] = 30

    if subresult_path == "geometric_graphs"
        title_string = "Random Geometric Graphs"
        save_title_string ="geometric_graphs"
    else
        title_string = "Erdős-Rényi Graphs"
        save_title_string ="erod_renyi_graphs"
    end

    if group_by == "density"
        bins = 1.0:-0.05:0.2
        y_label_string = "Graph Density"
    elseif group_by == "cutoff"
        bins = 1.0:-0.05:0.2
        y_label_string = "Cutoff"
    end

    for (i,graph_size) in enumerate(graph_sizes)
        df = collect_results(datadir("exp_raw/$subresult_path/$graph_size"))
        
        gdfs = groupby(df, :cutoff, sort=true)

        g_x_values = []
        g_y_values = []
        g_z_values = []
        for gdf in gdfs
            sys_nc_gdf = groupby(gdf, :max_inflight, sort=true)
            sys_df = sys_nc_gdf[1]
            nc_df = sys_nc_gdf[2]
            sort!(sys_df, [:graph_id])
            sort!(nc_df, [:graph_id])


            max_common_sample = min(nrow(sys_df),nrow(nc_df))
            gains = [sys_df.aoi_components[i][1]/nc_df.aoi_components[i][1] for i = 1:max_common_sample]
            
            if group_by == "density"
                append!(g_y_values,sys_df.cutoff[1:max_common_sample])
            elseif group_by == "cutoff"
                append!(g_y_values, rounds.([sys_df[i,:].cutoff for i in 1:max_common_sample],0.05))
            end
            append!(g_z_values, (gains.-1).*100)
        end
        gp_z_values = []
        gp_sample_values = []
        for bin in bins
            samples_indices = findall(g_y_values .== bin)
            if !isempty(samples_indices)
                append!(gp_z_values, mean(g_z_values[samples_indices]))
                append!(gp_sample_values, length(g_z_values[samples_indices]))
            else
                append!(gp_z_values, 0)
                append!(gp_sample_values, 0)
            end

        end

        push!(x_values, graph_size)
        push!(y_values, bins)
        push!(z_values, gp_z_values)
        push!(sample_values, gp_sample_values)
    end
    
    x_dim = length(x_values)
    y_dim = maximum([length(z) for z in z_values])
    
    data = zeros(y_dim,x_dim)
    labels = zeros(y_dim,x_dim)
    for x in 1:length(z_values)
        for y in 1:length(z_values[x])
            data[y,x] = z_values[x][y]
            labels[y,x] = sample_values[x][y]
        end
    end
    fig, ax = PyPlot.subplots()
    xlabel("\$N\$")
    ylabel(y_label_string)
    im = ax.imshow(data,cmap="YlOrBr")
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.94, anchor=(0,0), pad=0.06)
    cbar.ax.set_title("\$\\mathcal{G}_{nc}~[\\%]\$")

    ax.set_xticks((1:x_dim) .-1, labels=x_values)
    label_set = string.(collect(1.0:-0.05:0.2))
    label_set[2:2:end] .= ""
    ax.set_yticks((1:1:y_dim) .-1, labels=label_set)

    for i in 1:x_dim
        for j in 1:y_dim
            color = data[j,i]/maximum(data) > 0.75 ? [254,254,254]./255 : [0.0,0.0,0.0]
            color = "#$(hex(RGB(color[1],color[2],color[3])))"
            ax.text(i-1, j-1, "$(trunc(Int8,data[j, i]))", ha="center", va="center", color=color, fontsize=8, linespacing=1.2)
        end
    end
    im.set_clim(0,28)
    title("$title_string")
    savefig(plotsdir("graphsize_degree_gain-$save_title_string-$group_by.png"),bbox_inches="tight")
    savefig(plotsdir("graphsize_degree_gain-$save_title_string-$group_by.eps"),bbox_inches="tight")
end


pl_graphsize_degree_gain([8,10,12,14,16,18,20], "cutoff", "geometric_graphs")
pl_graphsize_degree_gain([8,10,12,14,16,18,20], "density", "erdos_renyi_graphs")