using
    Test,
    TWWL,
    CSV,
    DataFrames, DataFramesMeta,
    ProgressMeter,
    Logging, LoggingExtras

global_logger(ActiveFilteredLogger(log_args -> log_args.level != Logging.Warn, ConsoleLogger(stderr)))

@testset "Experiment1" begin
    testname = "Experiment1"
    if !(testname in ARGS) && length(ARGS) != 0
        return
    end

    nTrials = 10
    H = 5
    methods = ["Twwl", "WWL_1.0", "WWL_0.01", "LP"]
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
                trial=Int[],
                time=Float64[]
            )
        end

    for method in methods
        # for JIT. from here
        inp = downloadDataset("MUTAG")
        kernel =
            method == "Twwl" ? Twwl(H) :
            method == "LP" ? TwwlPrim(H) : DiscWWL(H)
        fit!(kernel, inp[:data])
        if method == "Twwl"
            treeMetric(kernel)
        elseif method == "LP"
            treeMetricLP(kernel)
        else
            wwlDist(kernel)
        end
        # to here

        for dataset in datasets
            if sum(@with df (:method .== method .&& :dataset .== dataset)) == nTrials
                continue
            end
            inp = downloadDataset(dataset)
            TWWL._updateUniqueMeasure!(inp[:data])

            kernels = Vector{Union{Twwl,TwwlPrim,DiscWWL}}(undef, nTrials)
            pts = zeros(Float64, nTrials)
            mts = zeros(Float64, nTrials)
            γs = zeros(Float64, nTrials)
            kernel =
                method == "Twwl" ? Twwl(H) :
                method == "LP" ? TwwlPrim(H) : DiscWWL(H)
            fit!(kernel, inp[:data]; verbose=true)
            for trial in 1:nTrials
                if any(@with df (:trial .== trial .&& :dataset .== dataset .&& :method .== method))
                    @info "[Skipped] trial: $(trial), dataset: $(dataset), method: $(method)"
                    continue
                end

                if method == "Twwl"
                    mts[trial] = TWWL.@timeout TIMEOUT (@elapsed treeMetric(kernel; verbose=true)) -1
                elseif method == "LP"
                    mts[trial] = TWWL.@timeout TIMEOUT (@elapsed treeMetricLP(kernel; verbose=true)) -1
                elseif startswith(method, "WWL")
                    γs[trial] =
                        method == "WWL_0.01" ? 0.01 :
                        method == "WWL_1.0" ? 1.0 : 10.0
                    mts[trial] = TWWL.@timeout TIMEOUT (@elapsed wwlDist(kernel; verbose=true, γ=γs[trial])) -1
                end
                push!(df, (dataset, method, trial, mts[trial]))
                CSV.write(outfile, df)
            end
        end
    end
end
