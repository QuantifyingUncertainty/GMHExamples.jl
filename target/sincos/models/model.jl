###Target functions
@everywhere function sincos!(r::Matrix,t::AbstractVector,paras::Vector)
    p1 = 2*pi*paras[1]
    p2 = 2*pi*paras[2]
    @simd for i=1:length(t)
        @inbounds r[i,1] = sin(p1*t[i])
    end
    @simd for i=1:length(t)
        @inbounds r[i,2] = cos(p2*t[i])
    end
    r
end

@everywhere sincos(t::AbstractVector,paras::Vector) = sincos!(zeros(eltype(t),length(t),2),t,paras)

@inline function _commonsincoscomponents(t::AbstractVector,measurementparas::Vector,variance::Vector,paraminit...)
    modelparameters = parameters([:a,:b],paraminit...)
    noisemodel = noise(:gaussian,variance)
    measurementdata = data(:array,t,applynoise!(noisemodel,sincos(t,measurementparas)))
    modelparameters,measurementdata,noisemodel
end

###Create a target function model
function sincosmodel(t::AbstractVector,measurementparas::Vector,variance::Vector,paraminit...)
    model(:target,_commonsincoscomponents(t,measurementparas,variance,paraminit...)...,sincos;name="Sine-Cosine Model")
end

###Create a target function model
function sincosmodel!(t::AbstractVector,measurementparas::Vector,variance::Vector,paraminit...)
    model(:target!,_commonsincoscomponents(t,measurementparas,variance,paraminit...)...,sincos!;name="Sine-Cosine In-Place Model")
end
