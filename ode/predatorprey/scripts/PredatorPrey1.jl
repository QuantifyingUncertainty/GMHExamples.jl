println("===================================================================")
println("Beginning GMHExamples.jl/ode/predatorprey/scripts/PredatorPrey1.jl")
println("===================================================================")

import GeneralizedMetropolisHastings #pre-include to avoid race condition when running in parallel
import GMHExamples

println("Number of parallel processes running: ",nprocs())
println("Use addprocs() in the REPL if you want to run on more processes")

@everywhere using GeneralizedMetropolisHastings
@everywhere using GMHExamples
using PyPlot

#Standard M-H for nproposals == 1
#Generalized M-H for nproposals > 1
nproposals = 200

#MCMC iteration specifications
nburnin = 500
niterations = 1000
ntunerperiod = 20

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
figure(string("PredatorPrey1-",nprocs())) ; clf()
subplot(221)
plot(dataindex(m),measurements(m)[:,1];label="Prey")
plot(dataindex(m),measurements(m)[:,2];label="Predator")
xlabel("Time")
ylabel("Amplitude")
title("Predator-Prey measurement data")
grid("on")
legend(loc="upper right",fancybox="true")

###Create a Metropolis Sampler with normal proposal density
s1 = sampler(:mh,:normal,1.0,[0.01,1.0,0.01,1.0,0.01,0.01])
s2 = sampler(:adaptive,0.01,numparas(m))
println("=============================")
println("Samplers defined successfully")
println("=============================")
show(s1)
show(s2)

###Create a tuner that scales the proposal density
t1 = tuner(:scale,ntunerperiod,0.5,:erf)
t2 = tuner(:monitor,ntunerperiod)
println("======+====================")
println("Tuners defined successfully")
println("=======+===================")
show(t1)
show(t2)

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
@time c1 = run!(r,m,s1,t1)
@time c2 = run!(r,m,s2,t2)
println("==========================")
println("Completed MCMC simulation")
println("=========================")

###Show the results of the simulations
show(c1)
show(c2)

nparas = numparas(m)
meanparamvals1 = mean(samples(c1),2)
stdparamvals1 = std(samples(c1),2)
meanparamvals2 = mean(samples(c2),2)
stdparamvals2 = std(samples(c2),2)

println("Results of the MH MCMC simulation:")
for i=1:nparas
    println(" parameter $(parameters(m)[i].key):  mean = ",meanparamvals1[i]," std = ",stdparamvals1[i])
end

println("Results of the Adaptive MCMC simulation:")
for i=1:nparas
    println(" parameter $(parameters(m)[i].key):  mean = ",meanparamvals2[i]," std = ",stdparamvals2[i])
end

println("================")
println("Plotting results")
println("================")

###Plot the loglikelihood values across samples
###After an initial few low values, this should remain relatively high
subplot(222)
plot(1:numsamples(c1),logposterior(c1,m))
title("Normal Metropolis-Hastings")
xlabel("Samples")
ylabel("Log-Posterior")

subplot(224)
plot(1:numsamples(c2),logposterior(c2,m))
title("Adaptive Metropolis-Hastings")
xlabel("Samples")
ylabel("Log-Posterior")

subplot(223)
modeldata1 = evaluate!(m,vec(meanparamvals1))
modeldata2 = evaluate!(m,vec(meanparamvals2))
plot(dataindex(m),measurements(m)[:,1];label="Prey")
plot(dataindex(m),measurements(m)[:,2];label="Predator")
plot(dataindex(m),modeldata1[:,1];label="MH Prey")
plot(dataindex(m),modeldata1[:,2];label="MH Predator")
plot(dataindex(m),modeldata2[:,1];label="AMH Prey")
plot(dataindex(m),modeldata2[:,2];label="AMH Predator")
xlabel("Time")
ylabel("Amplitude")
title("Predator-Prey model data")
grid("on")
legend(loc="upper right",fancybox="true")

###Plot the histograms of parameter values
figure(string("PredatorPrey-MH-",nprocs())) ; clf()
for i=1:nparas
    subplot(230 + i)
    h = PyPlot.plt[:hist](vec(getindex(samples(c1),i,:)),20)
    grid("on")
    title("Parameter $(parameters(m)[i].key)")
    xlabel("Values")
    ylabel("Samples")
end

figure(string("PredatorPrey-AMH-",nprocs())) ; clf()
for i=1:nparas
    subplot(230 + i)
    h = PyPlot.plt[:hist](vec(getindex(samples(c2),i,:)),20)
    grid("on")
    title("Parameter $(parameters(m)[i].key)")
    xlabel("Values")
    ylabel("Samples")
end
