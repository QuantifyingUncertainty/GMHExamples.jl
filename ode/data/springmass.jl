###Function to generate simulation data for the spring-mass dynamic system
### t      timepoints to evaluate the ODE at
### y0     initial conditions [x,dx/dt]
### paras  [M,K]
### Return: measurements
###  Column 1: position
###  Column 2: velocity
###Use examples/ode/test/plot_springmass.jl to plot data

###Function that can generate date for the spring-mass dynamic system
function springmassdata(t,y0,paras)
    a = sqrt(paras[1]/paras[2])
    y = [y0[1]*cos(a*t)+y0[2]/a*sin(a*t) -a*y0[1]*sin(a*t)+y0[2]*cos(a*t)]
end

springmassnoisy(n::GeneralizedMetropolisHastings.AbstractNoiseModel,t,y0,paras) = applynoise!(n,springmassdata(t,y0,paras))


