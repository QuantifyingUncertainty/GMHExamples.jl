## GeneralizedMetropolisHastings Examples
Example experiments for the Generalized Metropolis-Hastings algorithm [(Calderhead, 2014)](#refs). 

The elementary models, Julia scripts and IJulia Notebooks provided in this repository explain how to use the [GeneralizedMetropolisHastings.jl](https://github.com/QuantifyingUncertainty/GeneralizedMetropolisHastings.jl) library to develop new MCMC experiments. 

To run these experiments in Amazon Web Services, JuliaBox or on your local machine, see: http://quantifyinguncertainty.github.io

Top-level package contents:

1. **ode**: contains examples in which the model is an ordinary differential equation.
2. **function**: contains examples in which the model is a function which can be called directly.

Examples are organised in self-contained folders. They contain some or all of the following subfolders:

- **model**: model-specific Julia code
- **data**: files containing measurement data
- **notebooks**: documented examples for MCMC experiments, to be used in an IJulia notebook server
- **scripts**: documented examples for MCMC experiments that can be run from a Julia command-line (REPL) session
	
###<a name="refs"/>References
Calderhead B. (2014), A general construction for parallelizing Metropolis-Hastings algorithms, PNAS, Vol: 111, Pages: 17408-17413
