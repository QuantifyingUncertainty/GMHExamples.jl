## GMH-Examples.jl
Example experiments for the Generalized Metropolis-Hastings algorithm [(Calderhead, 2014)](#refs). 

The Julia scripts and IJulia Notebooks provided in this repository explain how to use the [GeneralizedMetropolisHastings.jl](https://github.com/QuantifyingUncertainty/GeneralizedMetropolisHastings.jl) library to develop new MCMC experiments.

To run these experiments in Amazon Web Services, JuliaBox or on your local machine, see: http://quantifyinguncertainty.github.io

Top-level package contents:

1. **ode**: contains examples in which the model is an ordinary differential equation.
2. **function**: contains examples in which the model is a function which can be called directly.
3. **photo**: contains code specific to the photo-receptor model

Each example directory is structured in the same way. It contains folders:

- **data**: measurement data for the physical/biological processes to infer parameters from
- **models**: types and functions that are model-specific, e.g., a spring-mass model or a Gaussian model
- **notebooks**: documented examples for MCMC experiments, to be used in an IJulia notebook server
- **scripts**: documented examples for MCMC experiments that can be run from a Julia command-line (REPL) session
	
	
The **test** directory contains tests of the package. You can run the tests from the Julia command line by going to the top-level folder of the GMH-Examples.jl module and executing:

```
julia> include("GMH-Examples.jl")
julia> include("test/runtests.jl")
```
These tests should complete successfully.

If you receive an error when executing the first of the two commands above, then it is likely that your installation misses the JLD package required to load and save julia data files. You can install this by executing:
```
julia> Pkg.add("JLD")
julia> Pkg.update()
```

###<a name="refs"/>References
Calderhead B. (2014), A general construction for parallelizing Metropolis-Hastings algorithms, PNAS, Vol: 111, Pages: 17408-17413
