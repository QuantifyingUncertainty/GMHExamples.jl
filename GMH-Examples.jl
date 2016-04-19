module GMHExamples

import Compat
import JLD
import Distributions
import GeneralizedMetropolisHastings

import Base: show

import GeneralizedMetropolisHastings:
    AbstractParameter,
    trait,_trait,traitvalue,traittype,
    _policy,
    parameter,parameters,_parameters,
    data,datavalues,dataindex,
    noise,noisemodel,
    model,_model,numvalues,measurements,
    evaluate!,
    applynoise!

export
    photonsequence, lightinducedcurrent, #from datafunctions.jl
    paramkeys,parampriors,fixedbumpvalues,fixedbumpshape, #from PhotoReceptorPolicy.jl
    photoreceptor,numtimesteps,numvilli,numphotons,microvilliwithphotons, #from PhotoReceptor.jl
    bump,bump!,lightinducedcurrent!,lightinducedcurrent #from lightinducedcurrent.jl

include("ode/data/fitzhughnagumo.jl")
include("ode/data/springmass.jl")
include("ode/models/fitzhughnagumo.jl")
include("ode/models/springmass.jl")
include("function/data/sincos.jl")
include("function/models/sincos.jl")
include("photo/data/datafunctions.jl")
include("photo/models/util.jl")
include("photo/models/PhotoReceptor.jl")
include("photo/models/lightinducedcurrent.jl")
include("photo/models/PhotoReceptorPolicy.jl")
include("photo/models/PhotoReceptorModel.jl")

end # module
