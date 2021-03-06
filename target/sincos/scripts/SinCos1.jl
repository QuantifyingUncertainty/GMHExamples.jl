println("=============================================================================")
println("Beginning execution of script GMHExamples.jl/target/sincos/scripts/SinCos1.jl")
println("=============================================================================")

import GeneralizedMetropolisHastings #pre-include to avoid race condition when running in parallel
import GMHExamples #pre-include to avoid race conditions

println("Number of parallel processes running: ",nprocs())
println("Use addprocs() in the REPL if you want to run on more processes")

@everywhere using GeneralizedMetropolisHastings
@everywhere using GMHExamples
using PyPlot


nproposals = 30

#MCMC iteration specifications
nburnin = 400
niterations = 1000
ntunerperiod = 20

#Time points to simulate the sine-cosine model
timepoints = linspace(0.0,10.0,500)

###Values of the model parameters (sine and cosine coefficients)
a = 1/2
b = 2/3
lows = [0.2,0.3]
highs = [0.8,1.0]

###The variance of the normal noise on the input data
variance = [0.09,0.09]

println("==========================================")
println("Simulation parameters defined successfully")
println("==========================================")

###Create a Spring-Mass model with measurement data and ODE function
m = sincosmodel(timepoints,[a,b],variance,lows,highs)

###Show the model
println("==========================")
println("Model defined successfully")
println("==========================")
show(m)

###Plot the measurement data (simmulated data + noise)
figure(string("SinCos1-",nprocs())) ; clf()
subplot(221)
plot(dataindex(m),measurements(m)[:,1];label="sine")
plot(dataindex(m),measurements(m)[:,2];label="cosine")
xlabel("Time")
ylabel("Amplitude")
title("Sine-Cosine data")
grid("on")
legend(loc="upper right",fancybox="true")

###Create a Metropolis sampler with a Normal proposal density
s = sampler(:mh,:normal,0.1,ones(2))
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
p = policy(:mh,nproposals;initialize=:prior)
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
println("=========================")
println("Completed MCMC simulation")
println("=========================")

###Show the result of the simulations
show(c)

nparas = numparas(m)
meanparamvals = mean(samples(c),2)
stdparamvals = std(samples(c),2)

println("Results of the MCMC simulation:")
println(" mean a: ",meanparamvals[1])
println(" mean b: ",meanparamvals[2])

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

###Plot the histograms of a,b values
for i=1:nparas
  subplot(222 + i)
    h = PyPlot.plt[:hist](vec(getindex(samples(c),i,:)),20)
  grid("on")
  title("Parameter $(parameters(m)[i].key)")
  xlabel("Values")
  ylabel("Samples")
end

###Finally, plot the average model results in the data window
subplot(221)
modeldata = evaluate!(m,vec(meanparamvals))
plot(dataindex(m),modeldata[:,1];label="model sine")
plot(dataindex(m),modeldata[:,2];label="model cosine")
legend(loc="upper right",fancybox="true")
