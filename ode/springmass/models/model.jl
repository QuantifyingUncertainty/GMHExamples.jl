###Function to generate data for a spring-mass dynamic system
### t: timepoints to evaluate
### y0: initial conditions
### paras: values of the parameters (K,M)
function springmassdata(t::AbstractVector,y0::AbstractVector,paras::AbstractVector)
    a = sqrt(paras[1]/paras[2])
    y = [y0[1]*cos(a*t)+y0[2]/a*sin(a*t) -a*y0[1]*sin(a*t)+y0[2]*cos(a*t)]
end

###Function to generate noisy data for a spring-mass dynamic system
### t: timepoints to evaluate
### y0: initial conditions
### paras: values of the parameters (K,M)
### noisegenerator: an AbstractNoiseModel used to apply measurement noise to the dynamic system
springmassnoisydata(t::AbstractVector,y0::AbstractVector,paras::AbstractVector,noisegenerator::AbstractNoiseModel) = applynoise!(noisegenerator,springmassdata(t,y0,paras))

###ODE for the spring-mass dynamic system
###The second-order differentional equation a = -K/M * x
###is written as a system of 2 first-order differential equations
###Function arguments are as expected by the Sundials ODE solver package
### t: timepoints to evaluate the ODE
### y: the values of the variables (x,v=dx/dt)
### ydot: the derivate vaues (v=dx/dt,a=d2x/dt2)
### paras: the equation parameters (K,M)
@everywhere function springmassode(t,y,ydot,paras)

    ydot[1] = y[2]
    ydot[2] = -paras[1]/paras[2]*y[1] #-K/M*X

    nothing
end

###Convenience function to create a spring-mass ODEModel
### t: timepoints to evaluate the ODE
### y0: initial conditions of the dynamic system
### measurementparas: the values for the parameters for which measurement data is generated (the "ground truth" of the experiment)
### variance: the variance of measurement noise
### paraminit: specification of the model parameters (K,M)
function springmassmodel(t::AbstractVector,y0::AbstractVector,measurementparas::AbstractVector,variance::AbstractVector,paraminit...)
    modelparameters = parameters([:K,:M],paraminit...)
    noisemodel = noise(:gaussian,variance)
    measurementdata = data(:function,t,springmassnoisydata,t,y0,measurementparas,noisemodel)
    model(:ode,modelparameters,measurementdata,noisemodel,springmassode,y0,2,[1,2];name = "Spring-Mass ODE")
end
