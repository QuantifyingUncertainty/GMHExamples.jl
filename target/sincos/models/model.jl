###Target functions
function sincos!(r::Matrix,paras::Vector,t::AbstractVector)
    p1 = 2pi*paras[1]
    p2 = 2pi*paras[2]
    @simd for i=1:length(t)
        @inbounds r[i,1] = sin(p1*t[i])
    end
    @simd for i=1:length(t)
        @inbounds r[i,2] = cos(p2*t[i])
    end
    r
end

sincos(paras::Vector,t::AbstractVector) = sincos!(zeros(eltype(t),length(t),2),paras,t)

@inline function _commonsincoscomponents(t::AbstractVector,measurementparas::Vector,variance::Vector,paraminit...)
    modelparameters = parameters([:a,:b],paraminit...)
    noisemodel = noise(:gaussian,variance)
    measurementdata = data(:array,t,applynoise!(noisemodel,sincos(measurementparas,t)))
    modelparameters,measurementdata,noisemodel
end

###Create a target function model
function sincosmodel(t::AbstractVector,measurementparas::Vector,variance::Vector,paraminit...)
    model(:target,_commonsincoscomponents(t,measurementparas,variance,paraminit...)...,sincos,t;name="Sine-Cosine Model")
end

###Create a target function model
function sincosmodel!(t::AbstractVector,measurementparas::Vector,variance::Vector,paraminit...)
    model(:target!,_commonsincoscomponents(t,measurementparas,variance,paraminit...)...,sincos!,t;name="Sine-Cosine In-Place Model")
end
