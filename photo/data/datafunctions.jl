##########
###Photons
##########

###Sequence of uniform random numbers
photonsequence(r::AbstractVector,len::Int,taxis =0.001) = GeneralizedMetropolisHastings.data(:array,taxis*(1:len),rand(r,len))

###Load an existing photon sequence
function photonsequence(filename::AbstractString,taxis =0.001f0)
    photons = trunc(Int,JLD.load(filename,"photons"))
    GeneralizedMetropolisHastings.data(:array,taxis*(1:length(photons)),photons)
end

##########
###Current
##########

###Load an existing light-induced current sequence
function lightinducedcurrent(filename::AbstractString,taxis =0.001f0)
    lic = JLD.load(filename,"lightinducedcurrent")
    GeneralizedMetropolisHastings.data(:array,taxis*(1:length(lic)),lic)
end

