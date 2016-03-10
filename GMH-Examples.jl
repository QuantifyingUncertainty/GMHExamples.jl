module GMHExamples

using Compat
using GeneralizedMetropolisHastings

import JLD
import Distributions

import Base: show
import GeneralizedMetropolisHastings: _trait,_policy,_parameters,_model,evaluate!

export
    AbstractPhotoReceptor,PhotoReceptorModel,
    fixedbumpshape,
    photonsequence,lightinducedcurrent,
    photoreceptor,numsteps,numvilli,numphotons,
    bump!,bump,macrocurrent!,macrocurrent,filterphotons!,filterphotons,
    lightinducedcurrent!,lightinducedcurrent

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
