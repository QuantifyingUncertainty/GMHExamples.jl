##########
###Photons
##########

###Sequence of uniform random numbers
photonsequence(r::AbstractVector{Int},len::Int) = rand(r,len)

###Load an existing photon sequence
photonsequence(filename::AbstractString) = trunc(Int,JLD.load(filename,"photons"))

nothing


