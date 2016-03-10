###ODE for the spring-mass dynamic system
###The second-order differentional equation
### a = -K/M * x
###is written as a system of first order differential equations
###Function arguments are as requested by the Sundials ODE solver package
### t = timepoints to evaluate the ODE
### y = the values of the variables (x,v=dx/dt)
### ydot = the derivate vaues (v=dx/dt,a=d2x/dt2)
### paras = the equation parameters (K,M)
function springmassode(t,y,ydot,paras)
  ydot[1] = y[2]
  ydot[2] = -paras[1]/paras[2]*y[1] #-K/M*x
end

###Helper function to create a spring-mass MCModel
### time to evaluate the ODE
### initial of the dynamic system
### modelvals the values for the model parameters
### paramvals the specifiers for the MCMC parameters
### noisevar the measurement noise variance
function springmassmodel(time,initial,modelvals,variance,paraminit...)
    p = parameters([:K,:M],paraminit...)
    n = noise(:gaussian,variance)
    d = data(:function,time,springmassnoisy,n,time,initial,modelvals)
    model(:ode,p,d,n,springmassode,initial,2,[1,2];name = "Spring-Mass")
end
