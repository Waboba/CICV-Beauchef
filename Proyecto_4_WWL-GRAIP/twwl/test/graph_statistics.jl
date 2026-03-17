using
    Test,
    TWWL,
    Graphs, MetaGraphsNext,
    Statistics, Printf

@testset "Output Basic Info" begin
    if !("Output Basic Info" in ARGS)
        return
    end

    datasets = [
        "MUTAG", "PTC_MR", "ENZYMES", "PROTEINS_full", "DD", "NCI1", "COLLAB",
        "REDDIT-BINARY", "REDDIT-MULTI-5K", "REDDIT-MULTI-12K", "DBLP_v1", "github_stargazers"
    ]

    for dataset in datasets
        println("----------[$(lpad(dataset, 18, ' '))]----------")
        inp = downloadDataset(dataset)
        N = length(inp[:target])

        for c in inp[:target] |> unique |> sort
            n = sum(inp[:target] .== c)
            print("$(n)/")
        end
        println("")

        n̄ = mean([nv(mg) for mg in inp[:data]])
        Ē = mean([ne(mg) for mg in inp[:data]])
        println("N\tn̄\tĒ")
        println(@sprintf("%d\t%.2f\t%.2f", N, n̄, Ē))

        println("-------------------------------------\n")
    end
end
