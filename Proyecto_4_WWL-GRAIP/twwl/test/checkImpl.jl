using
    Test,
    TWWL

@testset "ImplementationCheck" begin
    testname = "ImplementationCheck"
    if !(testname in ARGS) && length(ARGS) != 0
        return
    end

    H = 7
    datasets = ["MUTAG", "PTC_MR", "ENZYMES"]
    ϵ = 1e-10

    for dataset in datasets
        inp = downloadDataset(dataset)

        twwl = Twwl(H)
        twwlPrim = TwwlPrim(H)

        fit!(twwl, inp[:data]; verbose=true)
        fit!(twwlPrim, inp[:data]; verbose=true)

        dist = treeMetric(twwl; verbose=true)
        distPrim = treeMetric(twwlPrim; verbose=true)
        distLP = treeMetricLP(twwlPrim; verbose=true)

        @test isapprox.(dist, distPrim; atol=ϵ) |> all
        @test isapprox.(dist, distLP; atol=ϵ) |> all
    end
end

@testset "Check Impl of DiscWL" begin
    if !("Check Impl of DiscWL" in ARGS) && length(ARGS) != 0
        return
    end

    H = 7
    datasets = ["MUTAG", "PTC_MR", "ENZYMES", "DD", "PROTEINS_full", "NCI1"]
    ϵ = 1e-10

    for dataset in datasets
        inp = downloadDataset(dataset)

        twwlPrim = TwwlPrim(H)
        wwl = DiscWWL(H)

        fit!(twwlPrim, inp[:data])
        fit!(wwl, inp[:data])

        for (seq1, seq2) in zip(twwlPrim.labelSeq, wwl.seqs)
            @test isapprox.(seq1, seq2; atol=ϵ) |> all
        end
        for (measure1, measure2) in zip(twwlPrim.measures, wwl.measures)
            @test isapprox.(measure1, measure2; atol=ϵ) |> all
        end
    end
end
