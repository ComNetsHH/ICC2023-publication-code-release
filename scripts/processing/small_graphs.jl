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

include("$(srcdir())/utilities/topology_utilities.jl")
include("$(srcdir())/utilities/aoi_utilities.jl")


function pl_aoi_components(graph_sizes)
    
    gain_update = "#CCBB44"
    gain_encoding = "#EE6677"
    gain_full = "#4477AA"

    full_set_gains = []
    full_set_t_diff = []

    fig, ax = PyPlot.subplots(2,2)
    fig.tight_layout(rect=(0.045,0,1,0.95))
    fig.subplots_adjust(wspace=0.05, hspace=0.2)
    
    ax = [ax[1],ax[3],ax[2],ax[4]]
    plot_labeling = ["a)","b)","c)","d)"]

    ## Set up LaTeX fonts
    plt.rcParams["text.usetex"] = true
    plt.rcParams["font.family"] = "serif"
    plt.rcParams["font.serifx"] =  ["Computer Modern Roman"]
    plt.rcParams["font.size"] = 30
    plt.rc("hatch", linewidth=3.25)

    for (i,graph_size) in enumerate(graph_sizes)
        df = collect_results(datadir("exp_raw/$graph_size"))
        sys_nc_gdf = groupby(df, :max_inflight, sort=true)
        sys_df = sys_nc_gdf[1]
        nc_df = sys_nc_gdf[2]
        sort!(sys_df, [:graph_id])
        sort!(nc_df, [:graph_id])
        
        g_aoi_full = [sys[1] for sys in sys_df.aoi_components]./[nc[1] for nc in nc_df.aoi_components] .- 1
        g_aoi_full_capped = ([sys[1] for sys in sys_df.aoi_components]./[nc[1] for nc in nc_df.aoi_components] .- 1).*0.8
        g_aoi_update = ([sys[2] for sys in sys_df.aoi_components] - [nc[2] for nc in nc_df.aoi_components]) ./ [nc[1] for nc in nc_df.aoi_components]
        g_aoi_encoding = ([sys[3] for sys in sys_df.aoi_components] - [nc[3] for nc in nc_df.aoi_components]) ./ [nc[1] for nc in nc_df.aoi_components]
    
        sorter = sortperm(g_aoi_full)
        g_aoi_full = g_aoi_full[sorter]
        g_aoi_full_capped = g_aoi_full_capped[sorter]
        g_aoi_update = g_aoi_update[sorter]
        g_aoi_encoding = g_aoi_encoding[sorter]

        g_bar_positions = collect(1:length(g_aoi_full))
        if i > 2
            bar_width = 1.0
        else
            bar_width = 0.8
        end
        ratio = 0.66
        ax[i].bar(g_bar_positions, g_aoi_update.*100, color=gain_update, width=bar_width, edgecolor = gain_update)        
        
        if i > 2
            ratio = 1.0
            ax[i].bar(g_bar_positions.+(bar_width*(1-ratio))/2, g_aoi_full.*100 , color=gain_full, width=bar_width*ratio, edgecolor = gain_full)
            ax[i].bar(g_bar_positions, g_aoi_encoding.*100, color=gain_encoding, width=bar_width, edgecolor = gain_encoding)
        else
            ax[i].bar(g_bar_positions.+(bar_width*(1-ratio))/2, g_aoi_full.*100 , color=gain_full, width=bar_width*ratio, edgecolor = gain_full)
            ax[i].bar(g_bar_positions, g_aoi_encoding.*100, color=gain_encoding, width=bar_width, edgecolor = gain_encoding)
        end

        plt.sca(ax[i])
        xticks([])
        ylim([-12,45])

        xlabel("$(plot_labeling[i]) \$N=$graph_size\$",fontsize=12)
        if mod(i,2) == 1
            ylabel("Gain [%]")
        else
            plt.setp(ax[i], yticklabels=[])
        end
        ax[i].set_axisbelow(true)
        ax[i].grid(which="major",  color="#e7e7e7ff", linestyle="-", linewidth=0.5)
        yticks(collect(-10:10:45), fontsize=12)
        append!(full_set_gains,g_aoi_full)
        append!(full_set_t_diff,nc_df.t-nc_df.t_star)
    end
    
    lines = []
    append!(lines, ax[1].bar(3, -100 ,bottom=-100, color=gain_update))
    append!(lines, ax[1].bar(3, -100 ,bottom=-100, color=gain_encoding))
    append!(lines, ax[1].bar(3, -100 ,bottom=-100, color=gain_full))
    
    labels = ["\$\\mathcal{T}\$",
              "\$\\mathcal{C}\$",
              "\$\\mathcal{G}_{nc} = \\mathcal{T} + \\mathcal{C} \$",
              "\$U^{2}_{\\Delta}+C^{2}_{\\Delta}\$"]
    
    fig.legend(lines, labels, loc = (0.28, 0.925), ncol=4, edgecolor="black", fontsize=12)

    savefig(plotsdir("aoi_components.png"),bbox_inches="tight")
    savefig(plotsdir("aoi_components.eps"),bbox_inches="tight")
    savefig(plotsdir("aoi_components.pdf"),bbox_inches="tight")
    savefig(plotsdir("aoi_components.svg"),bbox_inches="tight")
    println("Average gain: $(mean(full_set_gains)*100) %")
    println("Maximum gain: $(maximum(full_set_gains)*100) %")
    println("In $(count(full_set_t_diff.>0)) out of $(size(full_set_t_diff)[1]) schedules, non t_star was optimal.")
end

graph_sizes = [4,5,6,7]
pl_aoi_components(graph_sizes)