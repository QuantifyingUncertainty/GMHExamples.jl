###Policies that determine the options in the model evaluation
immutable LatencyCalculation <: GeneralizedMetropolisHastings.AbstractPolicyTrait
    trait::Symbol
    LatencyCalculation(s::Symbol) = (@assert in(s,[:deterministic,:lognormal]) ; new(s))
end

immutable RefractoryCalculation <: GeneralizedMetropolisHastings.AbstractPolicyTrait
    trait::Symbol
    RefractoryCalculation(s::Symbol) = (@assert in(s,[:deterministic,:lognormal]) ; new(s))
end

immutable BumpShape <: GeneralizedMetropolisHastings.AbstractPolicyTrait
    trait::Symbol
    BumpShape(s::Symbol) = (@assert in(s,[:fixed,:sample]) ; new(s))
end

immutable FixedBumpShapeParameters{T<:AbstractFloat}
    amplitude::T
    shape::T
    scale::T
    cutoff::T
end

_traitname(t::LatencyCalculation) = :latency
_traitname(t::RefractoryCalculation) = :refractory
_traitname(t::BumpShape) = :bump

_trait(::Type{Val{:latency}},t::Symbol) = LatencyCalculation(t)
_trait(::Type{Val{:refractory}},t::Symbol) = RefractoryCalculation(t)
_trait(::Type{Val{:bump}},t::Symbol) = BumpShape(t)

paramkeys(t::GeneralizedMetropolisHastings.AbstractPolicyTrait) = tuple([symbol(_traitname(t),x) for x in _paramkeys(t)]...)
parampriors(t::GeneralizedMetropolisHastings.AbstractPolicyTrait) = _parampriors(t)

function _paramkeys(t::Union{LatencyCalculation,RefractoryCalculation})
    if t.trait == :deterministic
        return (:value,)
    elseif t.trait == :lognormal
        return (:location,:scale)
    else
        return ()
    end
end

function _paramkeys(t::BumpShape)
    if t.trait == :sample
        return (:amplitude,:shape,:scale)
    else
        return ()
    end
end

function _parampriors(t::Union{LatencyCalculation,RefractoryCalculation})
    if t.trait == :deterministic
        return (:uniform,)
    elseif t.trait == :lognormal
        return (:uniform,:uniform)
    else
        return ()
    end
end

function _parampriors(t::BumpShape)
    if t.trait == :sample
        return (:uniform,:lognormal,:lognormal)
    else
        return ()
    end
end

fixedbumpshape(;amplitude =4.0f0,shape =3.0f0,scale =2.5f0,cutoff =0.025f0) = FixedBumpShapeParameters(amplitude,shape,scale,cutoff)
fixedbumpvalues(f::FixedBumpShapeParameters) = (f.amplitude,f.shape,f.scale,f.cutoff)

immutable PhotoReceptorPolicy <: GeneralizedMetropolisHastings.AbstractPolicy
    latency::LatencyCalculation
    refractory::RefractoryCalculation
    bump::BumpShape
    fixedbump::FixedBumpShapeParameters
    PhotoReceptorPolicy(l,r,b,d) = new(l,r,b,d)
end

function _policy(::Type{Val{:photoreceptor}};latency =:lognormal,refractory =:lognormal,bump =:fixed,fixedbump =fixedbumpshape())
    PhotoReceptorPolicy(trait(:latency,latency),trait(:refractory,refractory),trait(:bump,bump),fixedbump)
end

function show(io::IO,b::FixedBumpShapeParameters)
    println(io,"  FixedBumpShapeParameters with values:")
    println(io,"    amplitude = ",b.amplitude)
    println(io,"    shape = ",b.shape)
    println(io,"    scale = ",b.scale)
    println(io,"    cutoff = ",b.cutoff)
    nothing
end

function show(io::IO,p::PhotoReceptorPolicy)
    println(io,"PhotoReceptorPolicy with traits:")
    println(io,"  latency = ",traitvalue(p.latency))
    println(io,"  refractory = ",traitvalue(p.refractory))
    println(io,"  bump = ",traitvalue(p.bump))
    show(io,p.fixedbump)
    nothing
end
