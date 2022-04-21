# HiQGA.jl

![CI status](https://github.com/GeoscienceAustralia/HiQGA.jl/workflows/CI/badge.svg)

This is a general purpose package for spatial statistical inference, geophysical forward modeling, Bayesian inference and inversion (both determinstic and probabilistic).

Readily usable geophysical forward operators are to do with AEM, CSEM and MT physics (references underneath), **for which the time domain AEM codes are fairly production-ready**. The current EM modeling is in 1D, but the inversion framework is dimensionally agnostic (e.g., you can regress images). Adding your own geophysical operators is easy, keep reading [down here](#developing-hiqga-or-modifying-it-for-your-own-special-forward-physics).

This package implements both the nested (2-layer) and vanilla trans-dimensional Gaussian process algorithm as published in 
- [*Bayesian inversion using nested trans-dimensional Gaussian processes*, A. Ray, Geophysical Journal International, **226(1)**, 2021](https://doi.org/10.1093/gji/ggab114).
- [*Bayesian geophysical inversion with trans-dimensional Gaussian process machine learning*, A. Ray and D. Myer, Geophysical Journal International **217(3)**, 2019](https://doi.org/10.1093/gji/ggz111).
- There is also a flavour of within-bounds Gauss-Newton/Occam's inversion implemented. For SkyTEM AEM, this is fully functional, but for other forward propagators you will have to provide a Jacobian (the linearization of the forward operator).

## Installation
To install, in a perfect world we'd use Julia's `Pkg` REPL by hitting `]` to enter `pkg>` mode. Then enter the following, at the `pkg>` prompt:
```
pkg> add https://github.com/GeoscienceAustralia/HiQGA.jl.git
```
## Usage
Examples of how to use the package can be found in the `examples` directory. Simply `cd` to the relevant example directory and `include` the .`jl` files in the order they are named. If using VSCode make sure to do *Julia: Change to this Directory* from the three dots menu on the top right. The Markov Chain Monte Carlo sampler is configured to support parallel tempering on multiple CPUs - some of the examples accomplish this with Julia's built-in multiprocessing, and others use MPI in order to support inversions on HPC clusters that don't work with Julia's default SSH-based multiprocessing. The MPI examples require [MPI.jl](https://github.com/JuliaParallel/MPI.jl) and [MPIClusterManagers.jl](https://github.com/JuliaParallel/MPIClusterManagers.jl/), which are not installed as dependencies for this package, so you will need to ensure they are installed and configured correctly to run these examples. Please note that MPIClusterManagers.jl has issues with Julia <1.4.2, so please ensure you are using an up-to-date Julia version. 

Some example scripts have as a dependency [Revise.jl](https://github.com/timholy/Revise.jl) as we're still actively [developing this package](https://pkgdocs.julialang.org/v1/getting-started/), so you may need to install [Revise](https://github.com/timholy/Revise.jl) if not already installed. All Julia users should be developing with [Revise](https://github.com/timholy/Revise.jl)! After installation, to run the examples, simply clone the package separately (or download as a ZIP), navigate to the `examples` folder and run the scripts in their numerical order.

## Updating the package 
```
pkg> update HiQGA
```

## Developing HiQGA or modifying it for your own special forward physics
After you have `]add`ed HiQGA you can simply do: 
```
pkg>dev HiQGA
```
If you haven't added it already, you can do:
```
pkg>dev https://github.com/GeoscienceAustralia/HiQGA.jl.git
```
It will download to your `JULIA_PKG_DEVDIR`. 

Another way is to simply clone or download this repository to your `JULIA_PKG_DEVDIR`, rename the cloned directory `HiQGA` removing the `.jl` bit and do
```
pkg>dev HiQGA
```
[Here's a gist](https://gist.github.com/a2ray/8c2c55c25fee6647501b403886bbe64d) on adding your own module if you want to modify the source code. Alternatively, if you only want to use the sampling methods in `HiQGA.transD_GP` without contributing to the source (boo! j/k) [here's another gist](https://gist.github.com/a2ray/92a8c14483c21dda6ddf56685b95fbb8) which is more appropriate. These gists were written originally for a package called `transD_GP` so you will have to modify `using transD_GP` to `using HiQGA.transD_GP`. Documentation is important and we're working on improving it before a full-release. 

## Development setup on NCI
The preferred development and usage environment for HiQGA is [Visual Studio Code](https://code.visualstudio.com/), which provides interactive execution of Julia code through the [VSCode Julia extension](https://code.visualstudio.com/docs/languages/julia). To install VSCode on the National Computational Infrastructure (NCI), you need to extract the VSCode rpm package using the steps in [this gist](https://gist.github.com/a2ray/701347f703b72abb630d2521b43c5f22), to a location where your account has write access.

It is also useful to use Revise.jl to ensure changes to the package are immediately reflected in a running Julia REPL (this is the reason that Revise is a dependency on some example scripts as noted above). More information on a workflow to use Revise during development can be found [here](https://gist.github.com/a2ray/e593751b24e45f8160ba8041fb811680).

### References for AEM and CSEM physics 

- [Blatter, D., Key, K., Ray, A., Foley, N., Tulaczyk, S., & Auken, E. (2018). Trans-dimensional Bayesian inversion of airborne transient EM data from Taylor Glacier, Antarctica. Geophysical Journal International, 214(3)](https://doi.org/10.1093/gji/ggy255)

- [Ray, A., & Key, K. (2012). Bayesian inversion of marine CSEM data with a trans-dimensional self parametrizing algorithm. Geophysical Journal International, 191(3), 1135-1151.](https://doi.org/10.1111/j.1365-246X.2012.05677.x)
