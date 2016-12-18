module GMHExamples

import GeneralizedMetropolisHastings:
    parameters,noise,data,model,
    applynoise!,
    AbstractNoiseModel

export
    sincosmodel,sincosmodel!,springmassmodel,fitzhughnagumomodel,predatorpreymodel

include("target/sincos/models/model.jl")
include("ode/springmass/models/model.jl")
include("ode/fitzhughnagumo/data/dataset1.jl")
include("ode/fitzhughnagumo/models/model.jl")
include("ode/predatorprey/data/dataset1.jl")
include("ode/predatorprey/models/model.jl")


end
