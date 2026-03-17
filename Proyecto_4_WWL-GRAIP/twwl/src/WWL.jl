using
    Graphs, MetaGraphsNext,
    Distances,
    OptimalTransport,
    Printf,
    ProgressMeter

export
    DiscWWL,
    ContWWL,
    wwlDist

mutable struct DiscWWL
    H::Int
    N::Union{Nothing,Int}
    labelSymb::Symbol
    measureSymb::Symbol
    seqs::Union{Nothing,Vector{AbstractMatrix}}
    measures::Union{Nothing,Vector{AbstractVector}}
end
function DiscWWL(H::Int; labelSymb::Symbol=:label, measureSymb::Symbol=:measure)
    _checkH(H)
    return DiscWWL(H, nothing, labelSymb, measureSymb, nothing, nothing)
end
function Base.show(io::IO, kernel::DiscWWL)
    text = @sprintf("""
        ----------[Discrete WWL struct]----------
        H: %d
        labelSymb: :%s
        measureSymb: :%s
        -----------------------------------------""",
        kernel.H, String(kernel.labelSymb), String(kernel.measureSymb)
    )
    print(io, text)
end

mutable struct ContWWL
    H::Int
    N::Union{Nothing,Int}
    D::Union{Nothing,Int}
    attrSymb::Symbol
    measureSymb::Symbol
    seqs::Union{Nothing,Vector{AbstractMatrix}}
    measures::Union{Nothing,Vector{AbstractVector}}
end
function ContWWL(H::Int; attrSymb::Symbol=:attr, measureSymb::Symbol=:measure)
    _checkH(H)
    return ContWWL(H, nothing, nothing, attrSymb, measureSymb, nothing, nothing)
end
function Base.show(io::IO, kernel::ContWWL)
    text = @sprintf("""
        ----------[Continuous WWL struct]----------
        H: %d
        attrSymb: :%s
        measureSymb: :%s
        -----------------------------------------""",
        kernel.H, String(kernel.attrSymb), String(kernel.measureSymb)
    )
    print(io, text)
end

function fit!(
    kernel::Union{DiscWWL,ContWWL}, mgs::Vector{MetaGraphsNext.MetaGraph};
    updateMeasure::Bool=true, verbose::Bool=false
)
    if updateMeasure
        _updateUniqueMeasure!(mgs; measureSymb=kernel.measureSymb)
    else
        _checkNodeInfo(mgs, kernel.measureSymb)
    end
    _checkNodeInfo(mgs, isa(kernel, ContWWL) ? kernel.attrSymb : kernel.labelSymb)
    D = isa(kernel, ContWWL) ? _getAttrDim(mgs, kernel.attrSymb) : 0

    N = length(mgs)
    seqs = Vector{AbstractMatrix}()
    for mg in mgs
        seq =
            if isa(kernel, ContWWL)
                zeros(Float64, nv(mg), D * (kernel.H + 1))
            else
                zeros(UInt, nv(mg), kernel.H + 1)
            end
        for uSymb in labels(mg)
            u = code_for(mg, uSymb)
            if isa(kernel, ContWWL)
                seq[u, 1:D] .= mg[uSymb][kernel.attrSymb]
            else
                labelHash = hash(mg[uSymb][kernel.labelSymb])
                seq[u, 1] = labelHash
            end
        end

        for h in 1:kernel.H
            startId = D * (h - 1) + 1
            endId = D * h

            for uSymb in labels(mg)
                u = code_for(mg, uSymb)

                if isa(kernel, ContWWL)
                    for vSymb in neighbor_labels(mg, uSymb)
                        v = code_for(mg, vSymb)
                        eInfo = (uSymb, vSymb) in keys(mg.edge_data) ?
                                mg.edge_data[uSymb, vSymb] :
                                mg.edge_data[vSymb, uSymb]
                        wgt = mg.weight_function(eInfo)

                        seq[u, (startId+D):(endId+D)] .+= (wgt .* seq[v, startId:endId]) ./ degree(mg, u)
                    end
                    if degree(mg, u) > 0
                        seq[u, (startId+D):(endId+D)] = (
                            seq[u, startId:endId] + seq[u, (startId+D):(endId+D)]
                        ) ./ 2.0
                    else
                        seq[u, (startId+D):(endId+D)] .= seq[u, startId:endId]
                    end
                else
                    labelHash = seq[u, h]
                    labelHashes = seq[code_for.(Ref(mg), neighbor_labels(mg, uSymb)), h] |> sort
                    nxt = hash(labelHash, hash(labelHashes))
                    seq[u, h+1] = nxt
                end
            end
        end
        push!(seqs, seq)
    end

    kernel.N = N
    if isa(kernel, ContWWL)
        kernel.D = D
    end
    kernel.seqs = seqs
    kernel.measures = [
        Float64[mg[label_for(mg, u)][kernel.measureSymb] for u in 1:nv(mg)]
        for mg in mgs
    ]

    return
end

"""
    wwlDist(kernel::DiscWWL)

    This function computes the Wasserstein Weisfeiler-Lehman distance matrix for labelled graphs,
    by using the Sinkhorn algorithm.
"""
function wwlDist(
    kernel::Union{DiscWWL,ContWWL}; verbose::Bool=false,
    γ::Float64=0.01, maxiter::Int=1000,
)
    (isnothing(kernel.N) || isnothing(kernel.seqs) || isnothing(kernel.measures)) &&
        throw(ErrorException("`kernel` must be applied `fit!`"))

    N = kernel.N
    dist = zeros(Float64, N, N)
    pbar = Progress(fld(N * (N - 1), 2); enabled=verbose, desc="wwlDist")
    for i in 1:N
        seqi = kernel.seqs[i]
        μ = kernel.measures[i]
        for j in i+1:N
            seqj = kernel.seqs[j]
            ν = kernel.measures[j]

            costs = pairwise(
                isa(kernel, ContWWL) ? Euclidean() : Hamming(),
                seqi, seqj; dims=1
            )
            opt = sinkhorn2(μ, ν, costs, γ; maxiter=maxiter)
            dist[i, j] = opt
            dist[j, i] = opt
            next!(pbar; showvalues=[("i", i), ("j", j)])
        end
        yield()
    end

    (dist + dist') ./ 2
end
