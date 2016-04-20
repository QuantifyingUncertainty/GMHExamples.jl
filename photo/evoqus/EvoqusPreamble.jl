#Make sure the package is installed
try
    Pkg.clone("git://github.com/QuantifyingUncertainty/GeneralizedMetropolisHastings.jl")
catch
    warn("Package already installed in previous run")
end

try
    Pkg.clone("git://github.com/QuantifyingUncertainty/GMHModels.jl")
catch
    warn("Package already installed in previous run")
end

import GeneralizedMetropolisHastings
import GMHModels

#Run the script
include("photo/evoqus/PhotoReceptor1.jl")
