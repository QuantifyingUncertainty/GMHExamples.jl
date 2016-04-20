#Run the script
println("=============================================================================")
println("Beginning execution of script GMH-Examples.jl/photo/evoqus/PhotoReceptor1.jl")
println("=============================================================================")

println("Number of parallel processes running: ",nprocs())

using GeneralizedMetropolisHastings
using GMHExamples

###Print a message indicating that the GMH package has loaded correctly
print_gmh_module_loaded()

println("================================")
println("Initialize Simulation Parameters")
println("================================")

#Standard M-H for nproposals == 1
#Generalized M-H for nproposals > 1
nproposals = 10000

#MCMC iteration specifications
nburnin = 100
niterations = 100
ntunerperiod = 10

###Values of the model
numvilli1 = 30000

#specify the values that determine the priors on the parameters
latencylocation = (2.0,3.5) #uniform distribution with (low,high) values
latencyscale = (0.2,0.7) #uniform distribution with (low,high) values
refractorylocation = (4.0,6.0) #uniform distribution with (low,high) values
refractoryscale = (1.5,2.5) #uniform distribution with (low,high) values
bumpamplitude = (3.0,5.0) #uniform distribution with (low,high) values
bumpshape = (log(3.0),0.3) #lognormal distribution with (location,scale) values; variable can vary roughly between 2.0 and 4.0, but becomes increasingly penalized outside
bumpscale = (log(2.5),0.3) #lognormal distribution with (location,scale) values; variable can vary roughly between 1.5 and 3.5, but becomes increasingly penalized outside

photons1 = photonsequence("photo/data/naturallight.jld")
current1 = lightinducedcurrent("photo/data/naturallight.jld")

modelpolicy1 = policy(:photoreceptor) #4-parameter model with stochastic lognormal latency and refractory parameters and fixed bump parameters
params1 = parameters(:photoreceptor,modelpolicy1,latencylocation,latencyscale,refractorylocation,refractoryscale)

####Variance for the noise model, estimated from previous runs
variance1 = [3600.0]

println("==========================================")
println("Simulation parameters defined successfully")
println("==========================================")

###Create a PhotoReceptor model
model1 = model(:photoreceptor,params1,photons1,current1,variance1,numvilli1,modelpolicy1)

###Show the model
println("==========================")
println("Model defined successfully")
println("==========================")
show(model1)

###Create a Metropolis sampler with a Normal proposal density
propcov = [0.1,0.01,0.1,0.01]
sampler1 = sampler(:mh,:normal,0.1,propcov)
println("============================")
println("Sampler defined successfully")
println("============================")
show(sampler1)

###Create a tuner that scales the proposal density
#tuner1 = tuner(:scale,ntunerperiod,0.5,:erf)
tuner1 = tuner(:scale,ntunerperiod,0.5,:erf)
println("==========================")
println("Tuner defined successfully")
println("==========================")
show(tuner1)

###Create a Generalized Metropolis-Hastings runner (which will default to Standard MH when nproposals=1)
runnerpolicy1 = policy(:mh,nproposals;initialize=:prior)
runner1 = runner(runnerpolicy1,niterations,nproposals;numburnin = nburnin)
println("===========================")
println("Runner defined successfully")
println("===========================")
show(runner1)

###Run the MCMC (can take quite a bit of time)
println("=======================")
println("Run the MCMC simulation")
println("=======================")
@time chain1 = run!(runner1,model1,sampler1,tuner1)
println("=========================")
println("Completed MCMC simulation")
println("=========================")

###Show the result of the simulations
show(chain1)

nparas = numparas(model1)
meanparamvals = mean(samples(chain1),2)
stdparamvals = std(samples(chain1),2)

println("Results of the MCMC simulation:")
for i=1:nparas
    println("mean $(parameters(model1)[i].key): $(meanparamvals[i])")
    println("std $(parameters(model1)[i].key): $(stdparamvals[i])")
end
