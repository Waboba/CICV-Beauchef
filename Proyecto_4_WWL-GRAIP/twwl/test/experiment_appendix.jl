using
    Test,
    TWWL,
    CSV,
    DataFrames, DataFramesMeta,
    ProgressMeter,
    Logging, LoggingExtras

@testset "ExperimentAppendix" begin
    testname = "ExperimentAppendix"
    if !(testname in ARGS) && length(ARGS) != 0
        return
    end

    nTrials = 10
    H = 7
    methods = ["Twwl", "WWL_1.0"]
    datasets = [
        "MUTAG", "PTC_MR", "ENZYMES", "PROTEINS_full", "DD", "NCI1",
        "COLLAB", "REDDIT-BINARY", "REDDIT-MULTI-5K", "REDDIT-MULTI-12K", "DBLP_v1", "github_stargazers",
    ]
    outDir = "out"
    outfile = joinpath(outDir, "$(testname).csv")
    TIMEOUT = 60 * 60 * 24

    df =
        if isfile(outfile)
            CSV.read(outfile, DataFrame; types=Dict("dataset" => String, "method" => String))
        else
            DataFrame(
                dataset=String[],
                method=String[],
                h=Int[],
                trial=Int[],
                time=Float64[]
            )
        end

    for method in methods
        # for JIT
        inp = downloadDataset("MUTAG")
        kernel = method == "Twwl" ? Twwl(H) : DiscWWL(H)
        fit!(kernel, inp[:data])
        if method == "Twwl"
            treeMetric(kernel)
        else
            wwlDist(kernel)
        end

        for dataset in datasets
            for h in 1:H
                if sum(@with df (:method .== method .&& :dataset .== dataset .&& :h .== h)) == nTrials
                    continue
                end

                inp = downloadDataset(dataset)
                TWWL._updateUniqueMeasure!(inp[:data])
                kernel = method == "Twwl" ? Twwl(h) : DiscWWL(h)
                fit!(kernel, inp[:data]; verbose=true)

                mts = zeros(Float64, nTrials)

                for trial in 1:nTrials
                    if any(@with df (:trial .== trial .&& :dataset .== dataset .&& :method .== method .&& :h .== h))
                        @info "[Skipped] trial: $(trial), dataset: $(dataset), method: $(method), h: $(h)"
                        continue
                    end

                    if method == "Twwl"
                        mts[trial] = @elapsed treeMetric(kernel; verbose=true)
                    else
                        mts[trial] = @elapsed wwlDist(kernel; verbose=true, γ=1.0)
                    end
                    push!(df, (dataset, method, h, trial, mts[trial]))
                    CSV.write(outfile, df)
                end
            end
        end
    end
end
