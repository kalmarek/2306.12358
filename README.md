# Replication for [2306.12358](https://arxiv.org/abs/2306.12358)
This repository contains notebooks and standalone code to replicate the computational results of _Kazhdan constants for Chevalley groups over the integers_ by Marek Kaluba and Dawid Kielak.

Before exploring the notebooks you need to clone the main repository:

```bash
git clone https://github.com/kalmarek/2306.12358.git
```

The replicating notebooks are located in `2306.12358/notebooks` subdirectory.

# Installation

You should install [julia](https://julialang.org/) from the [official repository](https://julialang.org/downloads/). Then while located in `2306/12358` directory run `julia` in terminal and execute the following commands in julias command-line (REPL) to instantiate the environment for computations.

```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

Instantiation should install (among others):

* [`JuMP`](https://jump.dev/) package for mathematical programming,
* [Splitting conic solver](https://github.com/cvxgrp/scs),
* [`COSMO`](https://github.com/oxfordcontrol/COSMO.jl) solver, and
* `IntervalArithmetic.jl` package from ValidatedNumerics (https://juliaintervals.github.io/).

The environment uses [`Groups.jl`](https://github.com/kalmarek/Groups.jl), [`StarAlgebras.jl`](https://github.com/kalmarek/StarAlgebras.jl/), [`SymbolicWedderburn.jl`](https://github.com/kalmarek/SymbolicWedderburn.jl/) and [`PropertyT.jl`](https://github.com/kalmarek/PropertyT.jl/) (unregistered) packages.

# Replication

## Notebooks

After instantiating the environment a jupyter server should be launched form `2306.12358` directory by issuing from julia command-line (REPL).

```julia
using Pkg
Pkg.activate(".")
using IJulia
notebook(dir=pwd())
```

> During the first run, the user may be asked for installation of Jupyter program (a server for running this notebook) within miniconda environment, which will happen automatically after confirmation. To execute the commands in the notebook, one needs to navigate to `notebooks` subdirectory and click either of the notebooks.

## Scripts

The following scripts are included in `scripts`  subdirectory:

* `SLnZ_AdjA2.jl` can be used to check that $\operatorname{Adj}_{A_2} - \lambda \Delta$ is positive in $\mathbb{R} \operatorname{SL}_{n}(\mathbb{Z})$. Note: this script, using the new language of gradings by root systems, reproves that $\operatorname{Adj}_n - \lambda \Delta \geqslant 0$ in $\mathbb{R} \operatorname{SL}_n(\mathbb{Z})$ (a result of Kaluba-Kielak-Nowak [KKN](https://arxiv.org/abs/1812.03456)).
* `Sp2nZ_AdjC2.jl` can be used to check that $\operatorname{Adj}_{C_2} - \lambda \Delta$ is positive in $\mathbb{R} \operatorname{Sp}_{2n}(\mathbb{Z})$.
* `G2_Adj.jl` can be used to check that $\operatorname{Adj}_{G_2} - \lambda \Delta$ is positive in $\mathbb{R} G_2$.

Additionally the two next scripts can be used to prove property T by computing sum of squares decomposition for $\Delta^2 - \lambda \Delta$:

* `Sp2nZ_has_T.jl` for group $\operatorname{Sp}_{2n}(\mathbb{Z})$.
* `G2_has_T.jl` for group $G_2$.

All of these scripts accept the command-line optional arguments:

```bash
  -N N                  the degree/genus/etc. parameter for a group
                        (type: Int64, default: 3)
  -R, --halfradius HALFRADIUS
                        the halfradius on which perform the sum of
                        squares decomposition (type: Int64, default: 2)
  -u, --upper_bound UPPER_BOUND
                        set upper bound for the optimization problem to
                        speed-up the convergence (type: Float64, default: Inf)
  -h, --help            show this help message and exit
```

Thus the following invocations can be used to reprove the computational results in the paper.

* Theorem 3.6:
  * `julia --project=@. scripts/SLnZ_AdjA2.jl -N 3 -R 2` 
  * `julia --project=@. scripts/SLnZ_AdjA2.jl -N 3 -R 3` (Note: the default solvers parameters in the script are not suitable for `-R 3`; in particular one has to increase `max_iters` considerably)
* Theorem 3.10:
  * `julia --project=@. scripts/Sp2nZ_AdjC2.jl -N 2 -R 3`
* Theorem 3.12:
  * `julia --project=@. scripts/SpnZ_has_T.jl -N 2 -R 2`
  * `julia --project=@. scripts/SpnZ_has_T.jl -N 2 -R 3`
* Theorem 3.15:
  * `julia --project=@. scripts/SpnZ_Level.jl -N 3 -R 2`
* Theorem 3.17:
  * `julia --project=@. scripts/G2_has_T.jl -R 2` (option `-N` is ignored)
* Theorem 3.18:
  * `julia --project=@. scripts/G2_Adj.jl -R 3` (option `-N` is ignored)

## Citing



If you find yourself using or studying code in this repository please cite

```
@misc{kaluba2023kazhdan,
      title={Kazhdan constants for Chevalley groups over the integers}, 
      author={Marek Kaluba and Dawid Kielak},
      year={2023},
      eprint={2306.12358},
      archivePrefix={arXiv},
      primaryClass={math.GR}
}
```

