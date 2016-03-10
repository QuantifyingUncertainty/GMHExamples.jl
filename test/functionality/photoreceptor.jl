using PyPlot

numvilli = 30000
latencyvals = [2.7137,0.4568]
refracvals = [log(150),log(7)]
bumpvals = [4.0,3.0,2.5]
paramvals = vcat(latencyvals,refracvals,bumpvals)

variance1 = [10.0] #variance for the noise model

#specify the values that determine the priors on the parameters
latencylocation = (2.0,3.5) #uniform distribution with (low,high) values
latencyscale = (0.2,0.7) #uniform distribution with (low,high) values
refractorylocation = (4.0,6.0) #uniform distribution with (low,high) values
refractoryscale = (1.5,2.5) #uniform distribution with (low,high) values
bumpamplitude = (3.0,5.0) #uniform distribution with (low,high) values
bumpshape = (log(3.0),0.3) #lognormal distribution with (location,scale) values; variable can vary roughly between 2.0 and 4.0, but becomes increasingly penalized outside
bumpscale = (log(2.5),0.3) #lognormal distribution with (location,scale) values; variable can vary roughly between 1.5 and 3.5, but becomes increasingly penalized outside

photons1 = photonsequence("photo/data/naturallight.jld")
zerocurrent1 = data(:array,dataindex(photons1),zeros(Float32,numvalues(photons1)))

policy1 = policy(:photoreceptor,bump=:sample) #stochastic lognormal latency and refractory parameters, sampled bumpshape parameters (with lognormal priors, see above)
params1 = parameters(:photoreceptor,policy1,latencylocation,latencyscale,refractorylocation,refractoryscale,bumpamplitude,bumpshape,bumpscale)

numevaluations = 100
nummodels = 10

meanresults1 = zeros(numvalues(photons1),nummodels)
varresults1 = zeros(numvalues(photons1),nummodels)
for i = 1:nummodels
    println("Creating a model with a new distribution of photons over microvilli")
    @time model1 = model(:photoreceptor,params1,photons1,zerocurrent1,variance1,numvilli,policy1)
    println("Evaluating the model $numevaluations times")
    @time results1 = evaluate!(model1,paramvals,numevaluations)
    meanresults1[:,i] = mean(results1,2)
    varresults1[:,i] = var(results1,2)
end
mmresults1 = mean(meanresults1,2)
mvresults1 = mean(varresults1,2)

figure("functionality/photoreceptor")
clf()
subplot(221)
tindex = dataindex(photons1)
plot(tindex,datavalues(photons1))
xlabel("Time")
ylabel("Number of Photons")
xlim(tindex[1],tindex[end])
title("Photon Sequence")
grid("on")

subplot(223)
plot(tindex,meanresults1)
plot(tindex,mmresults1,linewidth=4)
xlabel("Time")
ylabel("Current (nA)")
xlim(tindex[1],tindex[end])
title("Average light-induced current")
grid("on")

subplot(224)
plot(tindex,sqrt(varresults1))
plot(tindex,mean(sqrt(varresults1),2),linewidth=4)
xlabel("Time")
ylabel("STD of Current (nA)")
xlim(tindex[1],tindex[end])
title("Standard deviation of light-induced current")
grid("on")


modelcurrent1 = data(:array,dataindex(photons1),mmresults1)
variance1 = [mean(mvresults1)]

println("Variance obtained from model: $variance1")

println("Creating a model with experimentally-obtained macrocurrent and variance for the noise")
@time model2 = model(:photoreceptor,params1,photons1,modelcurrent1,variance1,numvilli,policy1)
println("Evaluating the model $numevaluations times with loglikelihood calculation")
loglikelihood1 = zeros(numevaluations)
@time for i=1:numevaluations
    l1 = evaluate!(model2,paramvals)
    loglikelihood1[i] = loglikelihood(model2,l1)
end
loglikelihood2 = zeros(numevaluations)
paramvals2 = vcat([2.6,0.4,log(145),log(6.5)],bumpvals)
@time for i=1:numevaluations
    l2 = evaluate!(model2,paramvals2)
    loglikelihood2[i] = loglikelihood(model2,l2)
end

subplot(222)
plot(1:numevaluations,loglikelihood1,label="ll1")
plot(1:numevaluations,loglikelihood2;label="ll2")
xlabel("Evaluation")
ylabel("Log-Likelihood")
xlim(1,numevaluations)
title("Log-Likelihood across evaluations")
grid("on")
legend(loc="upper right",fancybox="true")
