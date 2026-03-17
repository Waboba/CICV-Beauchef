using Distributed

#=
addprocs(
    [("localhost", 2)];
    tunnel=true,
    topology=:master_worker,
    enable_threaded_blas=true,
    exeflags="--project=~/Gits/TWWL"
)
=#

include("graph_statistics.jl")
include("checkImpl.jl")
include("experiment1.jl")
include("experiment2.jl")
include("experiment_appendix.jl")
