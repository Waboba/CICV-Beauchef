# README

This repository contains the experimental code for our paper "Tree Structure for the Categorical Wasserstein Weisfeiler-Lehman Graph Kernel".

Please note that any bugs or issues in this program are the responsibility of KeishiS.
If you encounter any problems please report them via GitHub issues and we will address them.

## Requirements

- Julia
- Python # we recommend to use [uv](https://github.com/astral-sh/uv).

## How to use

You can apply our algorithm to some datasets using the following commands:

```julia
> julia --project=.
(TWWL) pkg> instantiate
julia> using TWWL
julia> inp = downloadDataset("MUTAG")
julia> H = 5
julia> kernel = Twwl(H)
julia> fit!(kernel, inp[:data])
julia> dist = treeMetric(kernel)
```

Moreover, you can reproduce the results of our experiments using `Pkg.test` command.
Please note that unless you make some minor modifications for parallelization, the execution time can be extremely long.

```julia
# for TWWL, WWL
julia> using Pkg
julia> Pkg.test(; test_args=["Experiment1"]) # for Runtime comparison
julia> Pkg.test(; test_args=["Experiment2"]) # for Classification comparison
julia> Pkg.test(; test_args=["ExperimentAppendix"]) # for further runtime comparison
```

```python
# Before executing the following commands, please run Julia's `Pkg.test(; test_args=["Experiment2"])`.
> cd test/others
> uv sync
> source .venv/bin/activate
> METHOD=WL FILENAME=Experiment2_WL.csv python experiment2.py
> METHOD=SP FILENAME=Experiment2_SP.csv python experiment2.py
> METHOD=GL FILENAME=Experiment2_GL.csv python experiment2.py
> METHOD=WL-OA FILENAME=Experiment2_WL-OA.csv python experiment2.py
> METHOD=MLG FILENAME=Experiment2_MLG.csv python experiment2.py

# for R-WL
> deactivate
> cd GenWL && source .venv/bin/activate
> FILENAME=Experiment2_R-WL.csv python experiment2.py

# for WLLB
> deactivate
> cd WL-distance && source .venv/bin/activate
> FILENAME=Experiment2_WLLB.csv python experiment2.py
```

## Reference

```
@ARTICLE{Sando2025-zb,
title = "Tree Structure for the Categorical Wasserstein Weisfeiler-Lehman Graph Kernel",
author = "Sando, Keishi and Le, Tam and Hino, Hideitsu",
journal = "Transactions on Machine Learning Research",
month =  nov,
year =  2025
}
```
