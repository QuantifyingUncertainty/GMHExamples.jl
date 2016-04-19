import Base.Test: @test,@test_approx_eq,@test_approx_eq_eps,@test_throws

@everywhere import GeneralizedMetropolisHastings:
    ParameterUnivariate,ParameterDefault,
    trait,traitvalue,traittype,policy,
    parameters,
    data,datavalues,dataindex,numvalues,
    model,evaluate!,
    loglikelihood

@everywhere import GMHExamples:
    photonsequence,
    fixedbumpshape,
    photoreceptor,numtimesteps,numvilli,numphotons,
    bump!,bump,macrocurrent!,macrocurrent,filterphotons!,filterphotons,lightinducedcurrent!,lightinducedcurrent



