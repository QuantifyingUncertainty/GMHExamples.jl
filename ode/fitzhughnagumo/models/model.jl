#include measurment data
include ("../data/dataset1.jl")

###ODE for the FitzHugh-Nagumo model: http://www.scholarpedia.org/article/FitzHugh-Nagumo_model
### t = timepoints to evaluate the ODE at
### y = the model's state variables (Membrane Potential,Refractory Variable)
### ydot = the derivate values
### paras = the equation parameters (a,b,c)
@everywhere function fitzhughnagumoode(t,y,ydot,paras)

    a = paras[1]
    b = paras[2]
    c = paras[3]

    ydot[1] = c*(y[1]-(y[1]^3)/3+y[2])
    ydot[2] = -(y[1]-a+b*y[2])/c

    nothing
end

###Convenience function to create a FitzHugh-Nagumo ODEModel
### y0: initial conditions of the dynamic system
### variance: the variance of the noise model
### paraminit: the default and boundary values of the parameters
function fitzhughnagumomodel(y0::AbstractVector,variance::AbstractVector,paraminit...)
    modelparameters = parameters([:a,:b,:c],paraminit...)
    noisemodel = noise(:gaussian,variance)
    measurementdata = data(:array,fitzhughnagumodataset1()...)
    model(:ode,modelparameters,measurementdata,noisemodel,fitzhughnagumoode,y0,2,[1,2];name = "FitzHugh-Nagumo")
end
