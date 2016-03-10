###ODE for the FitzHugh-Nagumo model: http://www.scholarpedia.org/article/FitzHugh-Nagumo_model
### t = timepoints to evaluate the ODE at
### y = the model's state variables (Membrane Potential,Refractory Variable)
### ydot = the derivate values
### paras = the equation parameters (a,b,c)
function fitzhughnagumoode(t,y,ydot,paras)

  a = paras[1]
  b = paras[2]
  c = paras[3]

  ydot[1] = c*(y[1]-(y[1]^3)/3+y[2])
  ydot[2] = -(y[1]-a+b*y[2])/c
end

###Helper function to create a FitzHugh-Nagumo MCModel
### initial     initial conditions of the dynamic system
### variance    the variance of the noise model
### paraminit   the default and boundary values of the parameters
function fitzhughnagumomodel(initial,variance,paraminit...)
    p = parameters([:a,:b,:c],paraminit...)
    n = noise(:gaussian,variance)
    d = data(:array,fitzhughnagumodata()...)
    model(:ode,p,d,n,fitzhughnagumoode,initial,2,[1,2];name = "FitzHugh-Nagumo")
end
