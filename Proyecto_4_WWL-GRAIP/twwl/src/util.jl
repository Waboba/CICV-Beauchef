using
    Artifacts,
    SHA,
    DelimitedFiles,
    Graphs, MetaGraphsNext,
    Downloads,
    Serialization
using Pkg.Artifacts: bind_artifact!, create_artifact

export
    downloadDataset,
    paramGrid,
    @timeout

const supportedDatasets = Dict{String,Dict{Symbol,String}}(
    "MUTAG" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/MUTAG.zip",
        :sha256 => "c419bdc853c367d2d83da4973c45100954ae15e10f5ae2cddde6ca431f8207f6"
    ),
    "PTC_MR" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/PTC_MR.zip",
        :sha256 => "5699a6d9f1bc5b3d71495f09ef50de53fa3e6bb24ead1150da678500229f5237"
    ),
    "NCI1" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/NCI1.zip",
        :sha256 => "10e1458f3bd9224f14e6d7627e74dcfd13e48d376d73935e7bd2900590ef1d82"
    ),
    "DD" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/DD.zip",
        :sha256 => "d033d7aeb1a48c4b2b47cf0390e7fc9671de70d98c8df5b11e458d2ec5515cd2"
    ),
    "PROTEINS_full" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/PROTEINS_full.zip",
        :sha256 => "3b7782403ce98754df3330a67e9b2aff32e69520aa1245bf515c48cc0119c562"
    ),
    "ENZYMES" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/ENZYMES.zip",
        :sha256 => "13d832eb6ffa084192daf6e5750250028a18437ee692c38d29a10cd60e18aaf4"
    ),
    "SYNTHETIC" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/SYNTHETIC.zip",
        :sha256 => "3c0344e0cd6518d8b3f52bf45d152b1a9a007a523f5e91fc2bda929cc353c84d"
    ),
    "SYNTHETICnew" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/SYNTHETICnew.zip",
        :sha256 => "07e27d6ff1c25d036df5bf3593b1bf4676f51e656eb982e903dea4718634ed5e"
    ),
    "COIL-RAG" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/COIL-RAG.zip",
        :sha256 => "47f8be52a845c414299a4c9e6346acc074c5438fc12fe90b7d8da5c4613f667f"
    ),
    "Synthie" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/Synthie.zip",
        :sha256 => "c5bada5ffe42b4a901d50e75f10c3f969fb962f94acaab0265f46097496154d5"
    ),
    "IMDB-BINARY" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/IMDB-BINARY.zip",
        :sha256 => "b291ec8b26d85c70faa2ba0a2433e1f407ed2ef5d0fc072d36b9a95e49a1bb27"
    ),
    "IMDB-MULTI" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/IMDB-MULTI.zip",
        :sha256 => "a4a302149ebf4c76fa1f0fb108baff89fcbf9d35de306b18f27a8419b9a1a690"
    ),
    "BZR" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/BZR.zip",
        :sha256 => "584cb76f63cf5d0459cabfa254603b4e89823de0c03c96d60ab338f2d71a8b55"
    ),
    "COX2" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/COX2.zip",
        :sha256 => "aad0f932851a9ee308b3b65bb23d9eb94d4072de9593e988090ae837fa17c430"
    ),
    "BZR_MD" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/BZR_MD.zip",
        :sha256 => "f2ae267e19de998358bddeb20517768380ffb0594c594becbbd11aed538014d6"
    ),
    "COX2_MD" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/COX2_MD.zip",
        :sha256 => "c3a2849b0a6131dbc1a0f7213421e58c0ba2314db5dfe6787c0792674f631d9f"
    ),
    "REDDIT-BINARY" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/REDDIT-BINARY.zip",
        :sha256 => "982dc64ddade42a6365ba9c93f5080d0b9d315df0df492b57fd9fbb7e40dbe16"
    ),
    "REDDIT-MULTI-5K" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/REDDIT-MULTI-5K.zip",
        :sha256 => "68587783c33d54dd6107a1fc101594fd2b9d1af4e99cafe733636f383ae0c8ee"
    ),
    "REDDIT-MULTI-12K" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/REDDIT-MULTI-12K.zip",
        :sha256 => "dc133043edae3df088a4bfbd09595f105c0a6e2df533662b6ddc7d3268b5b49b"
    ),
    "reddit_threads" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/reddit_threads.zip",
        :sha256 => "b13c8a3984774bfd1a1e03f399939cea14cb94bc99afdeb6a3c5f1775bff9325"
    ),
    "twitch_egos" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/twitch_egos.zip",
        :sha256 => "b691bae434805f69c730c632e39ecde515c53fadfe1a9486741438bae9aa3160"
    ),
    "DBLP_v1" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/DBLP_v1.zip",
        :sha256 => "67d8a383e8920e9f9e6d8afd55df4104619252f2d08772ee08ad8a76e57b547a"
    ),
    "COLLAB" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/COLLAB.zip",
        :sha256 => "cb5b772bc4f3a4690d0601378a3fc993b092aeafd63e2b0f2330137dda1cdcc6"
    ),
    "github_stargazers" => Dict(
        :uri => "https://www.chrsmrrs.com/graphkerneldatasets/github_stargazers.zip",
        :sha256 => "9a2828ff7aca2ed9ff745edb3076e31b97b747f1e8f60ef193a3a5e8163f4680"
    )
)

macro timeout(seconds, expr, fail)
    quote
        tsk = @task $(esc(expr))
        schedule(tsk)
        Timer($(esc(seconds))) do timer
            istaskdone(tsk) || schedule(tsk, InterruptException(), error=true)
        end
        try
            fetch(tsk)
        catch _
            $(esc(fail))
        end
    end
end

"""
    downloadDataset(dataset::String; force::Bool=false)

    This function downloads datasets from [TUDatasets](https://chrsmrrs.github.io/datasets/).
    The return value is a dictionary containing the following keys:
        :data => Vector{MetaGraph}. Each MetaGraph object represents a graph structure and contains the node labels.
        :target => Vector{String}. Each element is the label of the corresponding graph.
"""
function downloadDataset(dataset::String; force::Bool=false)
    if !haskey(supportedDatasets, dataset)
        throw(ErrorException("""
            `$(dataset)` is not supported.
            The list of supported datasets is as follows: $(keys(supportedDatasets))"""))
    end
    info = supportedDatasets[dataset]

    artifactToml = joinpath(@__DIR__, "../Artifacts.toml")
    artifactHash = artifact_hash(dataset, artifactToml)

    if isnothing(artifactHash) || !artifact_exists(artifactHash) || force
        if isnothing(artifactHash) || !artifact_exists(artifactHash)
            artifactHash = create_artifact() do artifactDir
                zipfile = joinpath(artifactDir, basename(info[:uri]))
                Downloads.download(info[:uri], zipfile)
            end
            bind_artifact!(artifactToml, dataset, artifactHash; force=true)
        else
            Downloads.download(info[:uri], joinpath(@artifact_str(dataset), basename(info[:uri])))
        end
    end

    basepath = @artifact_str(dataset)
    zipfile = joinpath(basepath, basename(info[:uri]))
    binaryfile = joinpath(basepath, "data.jls")

    if isfile(binaryfile) && !force
        ret = open(binaryfile, "r") do io
            deserialize(io)
        end
        return ret
    end

    indicatorfile = "$(dataset)_graph_indicator.txt"
    adjacencyfile = "$(dataset)_A.txt"
    graphlabelsfile = "$(dataset)_graph_labels.txt"
    nodelabelsfile = "$(dataset)_node_labels.txt"
    edgelabelsfile = "$(dataset)_edge_labels.txt"
    nodeattrfile = "$(dataset)_node_attributes.txt"
    edgeattrfile = "$(dataset)_edge_attributes.txt"
    DUMMYLABEL = "_dummy_"

    # check hash value of zipfile
    hashBytes = open(zipfile, "r") do io
        sha256(io)
    end |> bytes2hex
    if hashBytes != info[:sha256]
        throw(ErrorException("The hash value of the downloaded zip file does not match the expected value."))
    end

    if !isfile(joinpath(basepath, indicatorfile)) || force
        # extract files
        run(pipeline(`unzip -o $(zipfile) -d $(basepath)`, stdout=devnull))
        for file in readdir(joinpath(basepath, dataset))
            filepath = joinpath(basepath, dataset, file)
            mv(filepath, joinpath(basepath, file); force=true)
        end
        rm(joinpath(basepath, dataset))
    end

    # indicatorfileのline numberがnode idであり，各lineの値がそのnode idが属するgraph id
    nid2gid = Dict{Int,Int}()
    graphs = Dict{Int,MetaGraph}()
    graphsInfo = Dict{Int,String}()
    eid2edge = Dict{Int,Tuple{Int,Int}}()
    N = readdlm(joinpath(basepath, graphlabelsfile), Int) |> vec |> length

    # mapping between node-id and graph-id
    for (nid, gid) in enumerate(readdlm(joinpath(basepath, indicatorfile), Int) |> vec)
        nid2gid[nid] = gid
        if !haskey(graphs, gid)
            graphs[gid] = MetaGraph(
                Graph(); label_type=Symbol,
                vertex_data_type=Dict{Symbol,Any},
                edge_data_type=Dict{Symbol,Any},
                weight_function=ed -> haskey(ed, :weight) ? ed[:weight] : 1.0,
                graph_data="Graph_" * lpad(string(gid), floor(Int, log10(N) + 1), '0')
            )
        end
    end

    # insert edge info
    for (eid, edge) in readdlm(joinpath(basepath, adjacencyfile), ',', Int) |> eachrow |> enumerate
        uNid, vNid = edge[1], edge[2]
        if nid2gid[uNid] != nid2gid[vNid]
            throw(ErrorException("something wrong 02"))
        end
        gid = nid2gid[uNid]
        uSymb, vSymb = Symbol(uNid), Symbol(vNid)
        eid2edge[eid] = (uNid, vNid)

        if !haskey(graphs[gid], uSymb)
            graphs[gid][uSymb] = Dict()
        end
        if !haskey(graphs[gid], vSymb)
            graphs[gid][vSymb] = Dict()
        end

        if !haskey(graphs[gid], uSymb, vSymb)
            add_edge!(graphs[gid], uSymb, vSymb, Dict())
        end
    end

    # insert node label/attributes info
    if isfile(joinpath(basepath, nodelabelsfile))
        for (nid, label) in enumerate(readdlm(joinpath(basepath, nodelabelsfile), String) |> vec)
            gid = nid2gid[nid]
            symb = Symbol(nid)
            if !haskey(graphs[gid], symb)
                graphs[gid][symb] = Dict(:label => label)
            else
                graphs[gid][symb][:label] = label
            end
        end
    else # In the absence of node label information, we assign the same label to all vertices.
        for (_, g) in graphs
            for uSymb in labels(g)
                deg = neighbor_labels(g, uSymb) |> length
                g[uSymb] = Dict(:label => "$(DUMMYLABEL)_$(deg)")
            end
        end
    end
    if isfile(joinpath(basepath, nodeattrfile))
        for (nid, attr) in readdlm(joinpath(basepath, nodeattrfile), ',', Float64) |> eachrow |> enumerate
            gid = nid2gid[nid]
            symb = Symbol(nid)
            if haskey(graphs[gid], symb)
                graphs[gid][symb][:attr] = attr
            else
                graphs[gid][symb] = Dict(:attr => attr)
            end
        end
    end

    # insert edge label/attributes info
    if isfile(joinpath(basepath, edgelabelsfile))
        for (eid, label) in readdlm(joinpath(basepath, edgelabelsfile), String) |> enumerate
            u, v = eid2edge[eid]
            symbU, symbV = Symbol(u), Symbol(v)
            gid = nid2gid[u]
            if haskey(graphs[gid], symbU, symbV)
                graphs[gid][symbU, symbV][:label] = label
            else
                throw(ErrorException("something wrong 04"))
            end
        end
    end
    if isfile(joinpath(basepath, edgeattrfile))
        for (eid, attr) in readdlm(joinpath(basepath, edgeattrfile), ',', Float64) |> eachrow |> enumerate
            u, v = eid2edge[eid]
            symbU, symbV = Symbol(u), Symbol(v)
            gid = nid2gid[u]
            if haskey(graphs[gid], symbU, symbV)
                if length(attr) == 1
                    graphs[gid][symbU, symbV][:weight] = attr[1]
                else
                    graphs[gid][symbU, symbV][:attr] = attr
                end
            else
                throw(ErrorException("something wrong 05"))
            end
        end
    end

    # graph labels
    for (gid, label) in readdlm(joinpath(basepath, graphlabelsfile), String) |> enumerate
        graphsInfo[gid] = label
    end

    if keys(graphs) != keys(graphsInfo)
        throw(ErrorException("something wrong 06"))
    end

    gkeys = keys(graphs) |> collect |> sort


    ret = Dict(
        :data => MetaGraph[graphs[k] for k in gkeys],
        :target => [graphsInfo[k] for k in gkeys]
    )
    open(binaryfile, "w") do io
        serialize(io, ret)
    end

    return ret
end

function _updateUniqueMeasure!(mgs::Vector{MetaGraphsNext.MetaGraph}; measureSymb::Symbol=:measure)
    for mg in mgs
        val = 1.0 / nv(mg)
        for uSymb in labels(mg)
            mg[uSymb][measureSymb] = val
        end
    end
end

function _checkNodeInfo(mgs::Vector{MetaGraphsNext.MetaGraph}, symb::Symbol)
    for mg in mgs
        for uSymb in labels(mg)
            if !haskey(mg[uSymb], symb)
                throw(ErrorException("$(symb) is necessary"))
            end
        end
    end
end

function _checkH(H::Int)
    H < 0 && throw(ErrorException("H must be non-negative"))
end

function _getAttrDim(mgs::Vector{MetaGraphsNext.MetaGraph}, attrSymb::Symbol)
    D = 0
    for mg in mgs
        for uSymb in labels(mg)
            if haskey(mg[uSymb], attrSymb)
                if D == 0
                    D = length(mg[uSymb][attrSymb])
                elseif D != length(mg[uSymb][attrSymb])
                    throw(ErrorException("The dimension of `attr` must be consistent accross all vertices"))
                end
            else
                throw(ErrorException("Vertex $(uSymb) does not have a $(attrSymb) attribute"))
            end
        end
    end

    D
end

function paramGrid(ps::Dict{Symbol,AbstractVector})::Vector{Dict}
    symbs = keys(ps)
    values = [ps[symb] for symb in symbs]

    [Dict(zip(symbs, comb)) for comb in Iterators.product(values...)] |> vec
end
