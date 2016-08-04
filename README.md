## GeneralizedMetropolisHastings Examples
Example experiments for the Generalized Metropolis-Hastings algorithm [(Calderhead, 2014)](#refs). 

The Julia scripts and IJulia Notebooks provided in this repository explain how to use the [GeneralizedMetropolisHastings.jl](https://github.com/QuantifyingUncertainty/GeneralizedMetropolisHastings.jl) library to develop new MCMC experiments. You will also need to install the [GMHModels.jl](https://github.com/QuantifyingUncertainty/GMHModels.jl) repository containing model-specific code.

To run these experiments in Amazon Web Services, JuliaBox or on your local machine, see: http://quantifyinguncertainty.github.io

Top-level package contents:

1. **ode**: contains examples in which the model is an ordinary differential equation.
2. **function**: contains examples in which the model is a function which can be called directly.
3. **photo**: contains code specific to the photo-receptor model

Example directories can contain folders:

- **notebooks**: documented examples for MCMC experiments, to be used in an IJulia notebook server
- **scripts**: documented examples for MCMC experiments that can be run from a Julia command-line (REPL) session
- **evoqus**: scripts to run on Evoqus
	
###<a name="refs"/>References
Calderhead B. (2014), A general construction for parallelizing Metropolis-Hastings algorithms, PNAS, Vol: 111, Pages: 17408-17413
