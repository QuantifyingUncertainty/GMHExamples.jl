module GMHExamples

using Compat
using GeneralizedMetropolisHastings

import JLD
import Distributions

export
  photonsequence,
  PhotoReceptor,numsteps,numvilli,distributephotons,
  bump!,bump,macrocurrent!,macrocurrent,lic

include("ode/data/fitzhughnagumo.jl")
include("ode/data/springmass.jl")
include("ode/models/fitzhughnagumo.jl")
include("ode/models/springmass.jl")
include("photo/data/photons.jl")
include("photo/models/util.jl")
include("photo/models/photoreceptor.jl")
include("photo/models/lightinducedcurrent.jl")


end # module
