using
    TWWL,
    Statistics, LinearAlgebra,
    MLJ, MLJBase, LIBSVM,
    ProgressMeter,
    DataFrames, DataFramesMeta,
    Arrow,
    Logging, LoggingExtras

params = Dict(
    :TWWL => Dict(
        :H => 1:7, :λ => 10.0 .^ range(-4, 1),
        :C => 10.0 .^ range(-3, 3)
    ),
    :WWL => Dict(
        :H => 1:7, :λ => 10.0 .^ range(-4, 1),
        :C => 10.0 .^ range(-3, 3), :γ => [0.01, 0.05, 0.1, 0.2, 0.5, 1.0, 10.0]
    )
)

cache = Dict()
cachey = Dict{String,AbstractVector}()
function precompute2(method, dataset)
    inp = downloadDataset(dataset)
    cachey[dataset] = inp[:target]

    for h in params[method][:H]
        kernel =
            method == :TWWL ? Twwl(h) :
            method == :WWL ? DiscWWL(h) :
            throw(ErrorException("undefined method"))
        TWWL.fit!(kernel, inp[:data])
        if method == :TWWL
            dist = treeMetric(kernel; verbose=true)
            cache[(method, dataset, h)] = dist
        elseif method == :WWL
            for γ in params[method][:γ]
                dist = wwlDist(kernel; verbose=true, γ=γ)
                cache[(method, dataset, h, γ)] = dist
            end
        else
            throw(ErrorException("undefined method"))
        end
    end
end

function eval2(method, dataset, trial, train, test; ntunecv::Int=10, seed::Int=42)
    y = cachey[dataset]

    # Parameter Tuning
    pbar = Progress(
        length(paramGrid(params[method]));
        desc="[$(method), $(dataset), $(trial)]"
    )
    scores = Dict{Dict,Float64}()
    for param in paramGrid(params[method])
        dist =
            method == :TWWL ? cache[(method, dataset, param[:H])] :
            method == :WWL ? cache[(method, dataset, param[:H], param[:γ])] :
            throw(ErrorException("undefined method"))
        dist_train = dist[train, train]
        y_train = y[train]
        accs = zeros(Float64, ntunecv)
        cv = StratifiedCV(; nfolds=ntunecv, shuffle=true, rng=seed)
        idxs = MLJBase.train_test_pairs(cv, 1:length(y_train), y_train) |> collect
        Threads.@threads for iter in eachindex(idxs)
            trainIdx = idxs[iter][1]
            validIdx = idxs[iter][2]
            K_tr = exp.(-param[:λ] .* dist_train[trainIdx, trainIdx])
            K_val = exp.(-param[:λ] .* dist_train[trainIdx, validIdx])
            y_tr = y_train[trainIdx]
            y_val = y_train[validIdx]

            model = svmtrain(Symmetric(K_tr), y_tr; kernel=Kernel.Precomputed, cost=param[:C])
            ypred, _ = svmpredict(model, K_val)
            accs[iter] = accuracy(y_val, ypred)
        end
        scores[param] = mean(accs)
        next!(pbar; showvalues=[
            ("total", pbar.n),
            ("counter", pbar.core.counter)
        ])
    end

    # Get bestparam and calculate accuracy with test dataset
    (bestparam, _) = reduce((a, b) -> a[2] > b[2] ? a : b, scores)
    bestdist =
        method == :TWWL ? cache[(method, dataset, bestparam[:H])] :
        method == :WWL ? cache[(method, dataset, bestparam[:H], bestparam[:γ])] :
        throw(ErrorException("undefined method"))

    K_train = exp.(-bestparam[:λ] .* bestdist[train, train])
    K_test = exp.(-bestparam[:λ] .* bestdist[train, test])
    y_train = y[train]
    y_test = y[test]

    bestmodel = svmtrain(K_train, y_train; kernel=Kernel.Precomputed, cost=bestparam[:C])
    ypred, _ = svmpredict(bestmodel, K_test)

    accuracy(y_test, ypred)
end

@testset "Experiment2" begin
    testname = "Experiment2"
    if !(testname in ARGS) && length(ARGS) != 0
        return
    end
    global_logger(ActiveFilteredLogger(log_args -> log_args.level != Logging.Warn, ConsoleLogger(stderr)))

    methods = [:TWWL, :WWL]
    datasets = [
        "MUTAG", "PTC_MR", "ENZYMES", "PROTEINS_full", "DD", "NCI1", "COLLAB",
        "REDDIT-BINARY", "REDDIT-MULTI-5K", "REDDIT-MULTI-12K", "DBLP_v1", "github_stargazers",
    ]

    seed = 42
    ncv = 10
    TIMEOUT = 60 * 60 * 24

    outDir = "out"
    outfile = joinpath(outDir, "$(testname).csv")
    idxfile = joinpath(outDir, "idx.arrow")
    if !isdir(outDir)
        mkdir(outDir)
    end

    df =
        isfile(outfile) ?
        CSV.read(outfile, DataFrame; types=Dict("method" => String, "dataset" => String)) :
        DataFrame(
            method=String[],
            dataset=String[],
            trial=Int[],
            accuracy=Float64[]
        )

    idxdf = DataFrame(
        dataset=String[],
        trial=Int[],
        trainIdx=Vector{Int}[],
        testIdx=Vector{Int}[]
    )

    for method in methods
        for dataset in datasets
            if sum(@with df (:method .== String(method) .&& :dataset .== dataset)) == ncv
                continue
            end

            inp = downloadDataset(dataset)
            N = length(inp[:target])
            precompute2(method, dataset)
            cv = StratifiedCV(; nfolds=ncv, shuffle=true, rng=seed)
            for (trial, (train, test)) in enumerate(MLJBase.train_test_pairs(cv, 1:N, inp[:target]))
                push!(idxdf, (dataset, trial, train, test))
                if any(@with df (:method .== String(method) .&& :dataset .== dataset .&& :trial .== trial))
                    continue
                end

                acc = @timeout TIMEOUT eval2(method, dataset, trial, train, test; seed=seed) -1
                push!(df, (
                    String(method), dataset, trial, acc
                ))
                CSV.write(outfile, df)

                acc == -1 && break
            end
        end
    end
    Arrow.write(idxfile, idxdf)
end
