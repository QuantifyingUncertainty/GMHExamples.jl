println("===================================================================")
println("Beginning GMH-Examples.jl/ode/predatorprey/scripts/PredatorPrey1.jl")
println("===================================================================")

import GeneralizedMetropolisHastings #pre-include to avoid race condition when running in parallel

println("Number of parallel processes running: ",nprocs())
println("Use addprocs() in the REPL if you want to run on more processes")

@everywhere using GeneralizedMetropolisHastings
using PyPlot

###Include the model functions
include("../models/model.jl")

#Standard M-H for nproposals == 1
#Generalized M-H for nproposals > 1
nproposals = 200

#MCMC iteration specifications
nburnin = 1000
niterations = 10000
ntunerperiod = 100

###Initial conditions for the ODE (prey and predator populations)
initial = [50.0,5.0]

###Default values of the parameters and prior boundaries
#defaults = [0.36,107.38,0.87,53.36,0.64,0.27]
defaults = [0.4,107.0,0.9,53.0,0.7,0.3]
lows = zeros(6)
highs = 150*ones(6)

###The variance of the noise on the input data
variance = sqrt(10.0)*ones(2)

println("==========================================")
println("Simulation parameters defined successfully")
println("==========================================")

###Create a FitzHugh-Nagumo model with measurement data, ODE function and parameters with default values and priors
m = predatorpreymodel(initial,variance,lows,highs,defaults)

###Show the model
println("==========================")
println("Model defined successfully")
println("==========================")
show(m)

###Plot the measurement data (simmulated data + noise)
figure("PredatorPrey1")
subplot(221)
plot(dataindex(m),measurements(m)[:,1];label="Prey")
plot(dataindex(m),measurements(m)[:,2];label="Predator")
xlabel("Time")
ylabel("Amplitude")
title("Predator-Prey measurement data")
grid("on")
legend(loc="upper right",fancybox="true")

###Create a Metropolis Sampler with normal proposal density
s = sampler(:mh,:normal,1.0,[0.005,0.08,0.012,0.03,0.006,0.003])
println("============================")
println("Sampler defined successfully")
println("============================")
show(s)

###Create a tuner that scales the proposal density
t = tuner(:scale,ntunerperiod,0.5,:erf)
println("==========================")
println("Tuner defined successfully")
println("==========================")
show(t)

###Create a Generalized Metropolis-Hastings runner (which will default to Standard MH when nproposals=1)
p = policy(:mh,nproposals;initialize=:default)
r = runner(p,niterations,nproposals;numburnin = nburnin)
println("===========================")
println("Runner defined successfully")
println("===========================")
show(r)

###Run the MCMC (can take quite a bit of time)
println("=======================")
println("Run the MCMC simulation")
println("=======================")
@time c = run!(r,m,s,t)
println("==========================")
println("Completed MCMC simulation")
println("=========================")

###Show the results of the simulations
show(c)

nparas = numparas(m)
meanparamvals = mean(samples(c),2)
stdparamvals = std(samples(c),2)

println("Results of the MCMC simulation:")
for i=1:nparas
    println(" parameter $(parameters(m)[i].key):  mean = ",meanparamvals[i]," std = ",stdparamvals[i])
end

println("================")
println("Plotting results")
println("================")

###Plot the loglikelihood values across samples
###After an initial few low values, this should remain relatively high
subplot(222)
plot(1:numsamples(c),logposterior(c,m))
title("Log-Posterior values across samples")
xlabel("Samples")
ylabel("Log-Posterior")

subplot(223)
modeldata = evaluate!(m,vec(meanparamvals))
plot(dataindex(m),measurements(m)[:,1];label="Prey")
plot(dataindex(m),measurements(m)[:,2];label="Predator")
plot(dataindex(m),modeldata[:,1];label="Model Prey")
plot(dataindex(m),modeldata[:,2];label="Model Predator")
xlabel("Time")
ylabel("Amplitude")
title("Predator-Prey model data")
grid("on")
legend(loc="upper right",fancybox="true")

###Plot the histograms of parameter values
figure("PredatorPrey2")
for i=1:nparas
  subplot(230 + i)
    h = PyPlot.plt[:hist](sub(samples(c),i,:)',20)
  grid("on")
  title("Parameter $(parameters(m)[i].key)")
  xlabel("Values")
  ylabel("Samples")
end










