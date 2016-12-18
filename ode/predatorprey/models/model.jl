###ODE for a predator-prey model
### t = timepoints to evaluate the ODE at
### y = the model's state variables (prey population size - predator population size)
### ydot = the derivate values
### paras = the equation parameters (r,K,s,a,u,v)
function predatorpreyode(t,y,ydot,paras)

    r = paras[1]
    k = paras[2]
    s = paras[3]
    a = paras[4]
    u = paras[5]
    v = paras[6]

    tmp = y[1]*y[2]/(a + y[1])

    ydot[1] = r*y[1]*(1-y[1]/k) - s*tmp
    ydot[2] = u*tmp - v*y[2]

    nothing
end

###Convenience function to create a Predator-Prey ODEModel
### y0: initial conditions of the dynamic system
### variance: the variance of the noise model
### paraminit: the default and boundary values of the parameters
function predatorpreymodel(y0::AbstractVector,variance::AbstractVector,paraminit...)
    modelparameters = parameters([:r,:k,:s,:a,:u,:v],paraminit...)
    noisemodel = noise(:gaussian,variance)
    measurementdata = data(:array,predatorpreydataset1()...)
    model(:ode,modelparameters,measurementdata,noisemodel,predatorpreyode,y0,2,[1,2];
          name="Predator-Prey",abstol=10e-9,reltol=10e-9)
end

