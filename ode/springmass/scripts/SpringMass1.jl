println("==================================================================================")
println("Beginning execution of script GMHExamples.jl/ode/springmass/scripts/SpringMass1.jl")
println("==================================================================================")

import GeneralizedMetropolisHastings #pre-include to avoid race condition when running in parallel
import GMHExamples

println("Number of parallel processes running: ",nprocs())
println("Use addprocs() in the REPL if you want to run on more processes")

@everywhere using GeneralizedMetropolisHastings
@everywhere using GMHExamples
using PyPlot

#Standard M-H for nproposals == 1
#Generalized M-H for nproposals > 1
nproposals = 10

#MCMC iteration specifications
nburnin = 200
niterations = 1000
ntunerperiod = 20

#Time points to simulate the spring-mass ODE
timepoints = 0.0:0.1:10.0

###Initial conditions for the spring-mass ODE (position and speed)
initialposition = -1.0 #in meters
initialvelocity = 1.0 #in meters/second

###Values of the model parameters (spring stiffness K and mass M)
K = 50.0 #in Newton/meter
M = 10.0 #in kg
lows = [K-K/5,M-M/5]
highs = [K+K/5,M+M/5]

###The variance of the normal noise on the input data
variance = [0.01,0.09]

println("==========================================")
println("Simulation parameters defined successfully")
println("==========================================")

###Create a Spring-Mass model with measurement data and ODE function
m = springmassmodel(timepoints,[initialposition,initialvelocity],[K,M],variance,lows,highs)

###Show the model
println("==========================")
println("Model defined successfully")
println("==========================")
show(m)

###Plot the measurement data (simmulated data + noise)
figure("SpringMass1")
subplot(221)
plot(dataindex(m),measurements(m)[:,1];label="location")
plot(dataindex(m),measurements(m)[:,2];label="velocity")
xlabel("Time")
ylabel("Amplitude")
title("Spring-Mass measurement data")
grid("on")
legend(loc="upper right",fancybox="true")

###Create a Metropolis sampler with a Normal proposal density
s = sampler(:mh,:normal,1.0,0.01eye(2))
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
println(" mean K: ",meanparamvals[1])
println(" mean M: ",meanparamvals[2])
println(" mean K/M: ",meanparamvals[1]/meanparamvals[2])
println("Mean K/M should be close to $(K/M)")

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

###Plot a scatter plot of K vs M values
###These should be spread around the K/M == 10.0 line (the diagonal in the figure)
ax3 = subplot(223)
ax3[:set_xlim]([lows[1],highs[1]])
ax3[:set_ylim]([lows[2],highs[2]])
scatter(vec(getindex(samples(c),1,:)),vec(getindex(samples(c),2,:)),marker=".",color="blue")
ax3[:set_aspect](abs(highs[1]-lows[1])/abs(highs[2]-lows[2]))
title("MCMC samples of Spring-Mass ODE parameters")
xlabel("Stiffness K (N/m)")
ylabel("Mass M (kg)")
grid("on")

###Plot a histogram of K/M values, which should peak around the true ratio of K/M
ax4 = subplot(224)
kml,kmu = K/M-K/M/20.0,K/M+K/M/20.0
ax4[:set_xlim]([kml,kmu])
nbins = linspace(kml,kmu,100)
h = PyPlot.plt[:hist]((vec(getindex(samples(c),1,:)./getindex(samples(c),2,:))),nbins)
grid("on")
xlabel("K/M")
ylabel("Number of Samples")
title("Histogram of K/M values")

###Finally, plot the average model results in the data window
subplot(221)
modeldata = evaluate!(m,vec(meanparamvals))
plot(dataindex(m),modeldata[:,1];label="model location")
plot(dataindex(m),modeldata[:,2];label="model velocity")






