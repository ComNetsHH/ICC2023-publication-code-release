[![DOI](https://zenodo.org/badge/ddd.svg)](https://zenodo.org/badge/latestdoi/ddd)

---

Any references to this code, protocol, or results shall be realized by citing the original paper from IEEE Xplore.

## Preface
This is a complementary code repository for the scientific publication 'Optimizing Age of Information in Status Update Systems using Network Coding: A Graph Search Approach' by Fisser, Leonard and Timm-Giel, Andreas published at the IEEE International Conference on Communications 2023 held in Rome, Italy.

---

## **Optimizing Age of Information in Status Update Systems using Network Coding — A Graph Search Approach**

Disseminating data in a wireless communication network is critical for the operation of future Cyber-Physical-Systems.
In the publication corresponding to this code repository, the problem of finding Age-of-Information optimal schedules for all-to-all data dissemination was investigated.
Age-optimal schedules were calculated by formulating a Graph-Search problem and deriving a custom A-Star search algorithm to solve for the dissemination schedule.
In addition to finding the already known optimal schedules for normal dissemination, the approach aims to find schedules which allow for Network Coding.
In Network Coding, multiple data packets may be combined on a bit-wise level to an encoded symbol carrying both information.
While Network Coding can increase the efficiency of a system, encoding and decoding delays may have adversary effects on the resulting Age-of-Information. 
As such, the proposed graph search approach balances the schedule reduction gains with coding delays and directly computes the age-optimal schedule given a specific topology.
Please refer to the original publication for more details.

## Getting Started
1. Download and install [JuliaLang](https://julialang.org/downloads/oldreleases/) Version 1.7.3.
2. Clone the repository to your local computer.
3. Start JuliaLang and navigate to the project's root directory.
4. Install project dependencies using:
```julia
pkg> activate .
pkg> instantiate
```
5. Familiarize yourself with the source code:
    - Calculations are started and configured by `scirpts/evaluation/analysis_runner.jl`
    - A-Star search variants are implemented in `src/solvers/a_star.jl` and `src/solvers/incremental_a_star.jl`
    - Plots and analysis are conducted using `DataFrames` and the code to reproduce the publication figures is located at `scripts/processing/`
    - Graph data and result data are stored in `data/exp_raw` and `data/exp_res` respectively. You may need to decompress the data sets used for the publication.
6. To run the calculations and find systematic / network coded schedules run the `analysis_runner.jl` via 
```julia
julia> include("scripts/evaluations/analysis_runner.jl")
```
7. You can access the raw result files using 
```julia
julia> df = collect_results(datadir("exp_raw/$graph_size"))
```

## Optional Packages
Both solvers need some general data and computations to set up the A-Star search algorithm.
Calculating this data may require the solving of linear programs.
These problems tend to get difficult for the inbuilt Julia solver `GLPK` with increasing graph sizes.
It is recommended to apply for a research license and use the [Gurobi](https://www.gurobi.com/) solver instead.
Search the project for mentions of `Gurobi` to find the code lines requiring a change.

## Comments and Notes
The provided code base was developed to answer specific research questions and its applicability and interfaces were not implemented with interoperability in mind.
Furthermore, the code was only tested on a single Linux machine and errors/bugs may occur on your system.
It is the intention of this repository to highlight how the results of the original publication were computed in detail.

## Special Thanks
The authors would like to thank the developers of the [DrWatson](https://github.com/JuliaDynamics/DrWatson.jl) JuliaLang package for their research project environment and helper functions. 

## Funding
This work was funded by the German Research Foundation (Deutsche Forschungsgemeinschaft, DFG) as part of the OUREL project with reference number 426655646.

