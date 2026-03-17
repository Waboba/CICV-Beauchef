module TWWL

using
    Printf,
    Graphs, MetaGraphsNext,
    SparseArrays,
    ProgressMeter,
    Artifacts

export
    Twwl,
    fit!,
    treeMetric

"""
    Twwl(H::Int)

    This returns a Twwl struct with the given maximum depth `H`.
    You can execute WL algorithm by `fit!` function.
"""
mutable struct Twwl
    H::Int
    N::Union{Nothing,Int}
    labelSymb::Symbol
    measureSymb::Symbol
    pIds::Union{Nothing,AbstractVector}
    measures::Union{Nothing,Dict{Int,AbstractVector}}
end
function Twwl(H::Int; labelSymb::Symbol=:label, measureSymb::Symbol=:measure)
    _checkH(H)
    return Twwl(H, nothing, labelSymb, measureSymb, nothing, nothing)
end
function Base.show(io::IO, kernel::Twwl)
    text = @sprintf("""
        ----------[Twwl struct]----------
        H: %d
        labelSymb: :%s
        measureSymb: :%s
        ---------------------------------""",
        kernel.H, String(kernel.labelSymb), String(kernel.measureSymb)
    )
    print(io, text)
end

function fit!(
    kernel::Twwl, mgs::Vector{MetaGraphsNext.MetaGraph};
    updateMeasure::Bool=true, verbose::Bool=false
)
    if updateMeasure
        _updateUniqueMeasure!(mgs; measureSymb=kernel.measureSymb)
    else
        _checkNodeInfo(mgs, kernel.measureSymb)
    end
    _checkNodeInfo(mgs, kernel.labelSymb)

    N = length(mgs)
    labelSeq = Dict{Int,AbstractMatrix}()
    pIds = sizehint!(Vector{Int}(), N * (kernel.H + 1))
    labelHash2nId = Dict{UInt,Int}()
    measures = Dict{Int,AbstractVector}()

    for (gId, mg) in enumerate(mgs)
        lseq = zeros(UInt, nv(mg), kernel.H + 1)
        for uSymb in labels(mg)
            u = code_for(mg, uSymb)
            labelHash = hash(mg[uSymb][kernel.labelSymb])
            lseq[u, 1] = labelHash

            if !haskey(labelHash2nId, labelHash)
                push!(pIds, 0)
                labelHash2nId[labelHash] = length(pIds)
            end
        end
        labelSeq[gId] = lseq
    end

    pbar = Progress((kernel.H + 1) * N; enabled=verbose, desc="fit!")
    for h in 0:kernel.H
        for (gId, mg) in enumerate(mgs)
            lseq = labelSeq[gId]
            for uSymb in labels(mg)
                u = code_for(mg, uSymb)
                labelHash = lseq[u, h+1]
                nId = labelHash2nId[labelHash]
                if h == kernel.H
                    if !haskey(measures, nId)
                        measures[nId] = spzeros(N)
                    end
                    measures[nId][gId] += mg[uSymb][kernel.measureSymb]
                else
                    labelHashes = lseq[code_for.(Ref(mg), neighbor_labels(mg, uSymb)), h+1] |> sort
                    nxt = hash(labelHash, hash(labelHashes))
                    lseq[u, h+2] = nxt

                    if !haskey(labelHash2nId, nxt)
                        push!(pIds, nId)
                        labelHash2nId[nxt] = length(pIds)
                    end
                end
            end
            next!(pbar)
        end
    end
    kernel.N = N
    kernel.pIds = pIds
    kernel.measures = measures

    return
end

"""
    treeMetric(kernel::Twwl)

    This function returns the Wasserstein Weisfeiler-Lehman distance matrix for labelled graphs.
    This function is based on the TWWL algorithm and is implemented in Julia.
"""
function treeMetric(kernel::Twwl; verbose::Bool=false)::AbstractMatrix
    (isnothing(kernel.pIds) || isnothing(kernel.measures) || isnothing(kernel.N)) &&
        throw(ErrorException("`kernel` must be applied `fit!`"))

    N = kernel.N
    dist = zeros(Float64, N, N)

    measures = deepcopy(kernel.measures)
    pbar = Progress(length(kernel.pIds); enabled=verbose, desc="treeMetric")
    for (nId, pId) in Iterators.reverse(enumerate(kernel.pIds))
        measure = pop!(measures, nId)
        nzIdx, vals = findnz(measure)
        dist[nzIdx, :] .+= vals ./ 2.0
        dist[nzIdx, nzIdx] .-= vals ./ 2.0

        dist[:, nzIdx] .+= vals' ./ 2.0
        dist[nzIdx, nzIdx] .-= vals' ./ 2.0

        dist[nzIdx, nzIdx] .+= abs.(vals .- vals') ./ 2.0

        if !haskey(measures, pId)
            measures[pId] = spzeros(N)
        end
        measures[pId] .+= measure
        next!(pbar)

        if nId % 100 == 0
            yield() # this is for cooperative multitasking
        end
    end

    (dist + dist') ./ 2
end

include("util.jl")
include("TwwlPrim.jl")
include("WWL.jl")

end
