println("============================================================================")
println("Beginning execution of script GMH-Examples.jl/ode/scripts/FitzHughNagumo1.jl")
println("============================================================================")

import GeneralizedMetropolisHastings
import GMHModels

println("Number of parallel processes running: ",nprocs())
println("Use addprocs() in the REPL if you want to run on more processes")

@everywhere using GeneralizedMetropolisHastings
@everywhere using GMHModels
using PyPlot

###Print a message indicating that the GMH package has loaded correctly
print_gmh_module_loaded()

println("================================")
println("Initialize Simulation Parameters")
println("================================")

#Standard M-H for nproposals == 1
#Generalized M-H for nproposals > 1
nproposals = 30

#MCMC iteration specifications
nburnin = div(3000,nproposals)
niterations = div(30000,nproposals)
ntunerperiod = div(300,nproposals)

###Initial conditions for the spring-mass ODE (membrane potential and refractory variable)
initial = [-1.0,1.0]

###Default values of the parameters (a,b,c) and prior boundaries
defaults = [0.3,0.3,2.0]
lows = zeros(3)
highs = [5.0,5.0,5.0]

###The variance of the noise on the input data
variance = [0.02,0.005]

println("==========================================")
println("Simulation parameters defined successfully")
println("==========================================")

###Create a FitzHugh-Nagumo model with measurement data, ODE function and parameters with default values and priors
m = fitzhughnagumomodel(initial,variance,lows,highs,defaults)

###Show the model
println("==========================")
println("Model defined successfully")
println("==========================")
show(m)

###Plot the measurement data (simmulated data + noise)
figure("FitzHughNagumo1")
subplot(231)
plot(dataindex(m),measurements(m)[:,1];label="membrane potential")
plot(dataindex(m),measurements(m)[:,2];label="refractory variable")
xlabel("Time")
ylabel("Amplitude")
title("FitzHugh-Nagumo measurement data")
grid("on")
legend(loc="upper right",fancybox="true")

###Create a Metropolis sampler with a Normal proposal density
s = sampler(:mh,:normal,0.01,eye(3))
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
subplot(232)
plot(1:numsamples(c),logposterior(c,m))
title("Log-Posterior values across samples")
xlabel("Samples")
ylabel("Log-Posterior")

###Plot the histograms of a,b,c values
for i=1:nparas
  subplot(233 + i)
    h = PyPlot.plt[:hist](sub(samples(c),i,:)',20)
  grid("on")
  title("Parameter $(parameters(m)[i].key)")
  xlabel("Values")
  ylabel("Samples")
end

subplot(233)
modeldata = evaluate!(m,vec(meanparamvals))
plot(dataindex(m),measurements(m)[:,1];label="membrane potential")
plot(dataindex(m),measurements(m)[:,2];label="refractory variable")
plot(dataindex(m),modeldata[:,1])
plot(dataindex(m),modeldata[:,2])
xlabel("Time")
ylabel("Amplitude")
title("FitzHugh-Nagumo model data")
grid("on")
legend(loc="upper right",fancybox="true")








