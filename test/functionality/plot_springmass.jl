### Load the plot packages
using Gadfly
using PyPlot

### Print a help message
println("Plot spring-mass measurement data using:")
println(" Gadfly package: gadflyplot_springmass(y0,paras,noisevar;timepoints)")
println(" PyPlot package: pyplot_springmass(y0,paras,noisevar;timepoints)")

### Plot the spring-mass data using the Gadfly package
function gadflyplot_springmass(timepoints,y0,paras,noisevar)

    #get the data
    measurements = GMHExamples.springmassnoisy(GeneralizedMetropolisHastings.noise(:gaussian,noisevar),timepoints,y0,paras)

    #plot the data
    Gadfly.plot(x=timepoints,y=measurements)

end

### Plot the spring-mass data using the PyPlot package
function pyplot_springmass(timepoints,y0,paras,noisevar)

    #get the data
    measurements = GMHExamples.springmassnoisy(GeneralizedMetropolisHastings.noise(:gaussian,noisevar),timepoints,y0,paras)

    #plot the data
    fig = PyPlot.figure()
    PyPlot.plot(timepoints,measurements[:,1];label="Position")
    PyPlot.plot(timepoints,measurements[:,2];label="Velocity")
    PyPlot.xlabel("Time")
    PyPlot.ylabel("Amplitude")
    PyPlot.title("Spring-Mass Measurement Data")
    PyPlot.grid("on")
    PyPlot.legend(loc="upper right",fancybox="true")

    #return the figure object
    fig
end

pyplot_springmass(0.0:0.1:10.0,[1,-1],[50.0,10.0],[0.1,0.2])

# File return statement
nothing
