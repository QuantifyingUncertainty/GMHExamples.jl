nsteps = 5
nvilli = 8
rphotons = 0:10

###Test creation of photonsequence
photons = photonsequence(rphotons,nsteps)
nphotons = sum(datavalues(photons))
@test numvalues(photons) == nsteps
@test all(minimum(rphotons) .<= datavalues(photons) .<= maximum(rphotons))

###Test utility functions
@test GMHExamples.uniquesorted(Int[]) == Int[]
@test GMHExamples.uniquesorted([2]) == [2]
@test GMHExamples.uniquesorted([1,1,1]) == [1]
@test GMHExamples.uniquesorted([1,2,2,3,3,3,8]) == [1,2,3,8]
@test GMHExamples.uniquesorted([1,4,7]) == [1,4,7]

@test_approx_eq (srand(465) ; GMHExamples.randlognormal!(zeros(5),1.0,2.0)) (srand(465) ; [rand(Distributions.LogNormal(1.0,2.0)) for i=1:5])
@test_approx_eq (srand(465) ; GMHExamples.randlognormal!(zeros(Float32,5),1.0f0,2.0f0)) (srand(465) ; [Float32(rand(Distributions.LogNormal(1.0,2.0))) for i=1:5])

###Test microvilli implementations
@test_throws MethodError photoreceptor(:dummy,datavalues(photons),nvilli)
receptors = Dict{Symbol,Any}()
for args in [(:full,(Int,)),(:row,(UInt8,Int16)),(:cell,(UInt16,))]
    srand(987)
    receptors[args[1]] = photoreceptor(args[1],datavalues(photons),nvilli,args[2]...)
    @test numvilli(receptors[args[1]]) == nvilli
    @test numsteps(receptors[args[1]]) == nsteps
    @test numphotons(receptors[args[1]]) == nphotons
end

@test typeof(receptors[:full]) <: GMHExamples.PhotoReceptorMatrix
@test typeof(receptors[:row]) <: GMHExamples.PhotoReceptorMatrix
@test typeof(receptors[:cell]) <: GMHExamples.PhotoReceptorCell

for t=1:nsteps
  @test GMHExamples.microvilliwithphotons(receptors[:full],t) == GMHExamples.microvilliwithphotons(receptors[:row],t)
  @test GMHExamples.microvilliwithphotons(receptors[:full],t) == unique(GMHExamples.microvilliwithphotons(receptors[:cell],t))
end

###Test current integration functions
@test GMHExamples._bumpvalmaximum(0.0,0.0) == 0.0
@test GMHExamples._bumpvalmaximum(3.0,2.5) == 7.5
@test GMHExamples._bumpvalmaximum(3.0f0,2.5f0) == 7.5f0 #Float32
@test GMHExamples._bumpval(1.0,1.0,1.0,0.0) == 0.0
@test GMHExamples._bumpval(2.0,3.0,1.5,0.0) == 0.0
@test GMHExamples._bumpval(2.0,3.0,1.5,4.5) == 2.0
@test GMHExamples._bumpval(4.0f0,3.0f0,2.5f0,7.5f0) == 4.0f0 #Float32
@test bump!(Vector{Float32}(30),collect(1.0f0:1.0f0:30.0f0),4.0f0,3.0f0,2.5f0) == bump(collect(1.0f0:1.0f0:30.0f0),4.0f0,3.0f0,2.5f0)
@test bump!(Vector{Float64}(30),collect(1.0:1.0:30.0),4.0,3.0,2.5) == bump(collect(1.0:1.0:30.0),4.0,3.0,2.5)
@test length(bump(0.0,0.0,0.0,0.01)) == 0

bmp1 = bump(1.0,1.0,1.0,0.025)
@test maximum(bmp1) == 1.0 && GMHExamples._bumpvalmaximum(1.0,1.0) == 1.0 && length(bmp1) == 6

ph = zeros(Int,8,7) ; mc = zeros(8,7)
ph[1,1] = 1 ; mc[1:6,1] = bmp1
ph[2,2] = 1 ; mc[2:7,2] = bmp1
ph[4,3] = 1 ; mc[4:8,3] = bmp1[1:5]
ph[1,4] = 2 ; mc[:,4] = 2*mc[:,1]
ph[1:2,5] = [1,1] ; mc[:,5] = mc[:,1] + mc[:,2]
ph[[2,4],6] = [1,1] ; mc[:,6] = mc[:,2] + mc[:,3]
ph[[1,2,4],7] = [2,1,1] ; mc[:,7] = mc[:,4] + mc[:,6]

for j=1:size(ph,2)
  @test_approx_eq mc[:,j] macrocurrent!(zeros(8),ph[:,j],bmp1)
  @test_approx_eq mc[:,j] macrocurrent(ph[:,j],bmp1)
end
@test_throws AssertionError macrocurrent!(zeros(10),ph[:,1],bmp1)

receptor1 = GMHExamples.PhotoReceptorMatrix{Int,Matrix}([1 0 1 0 0 1 1 1;1 0 0 1 1 0 0 1;0 0 0 0 0 0 0 1]',8,3,10)

@test filterphotons(receptor1,ones(sum(receptor1.microvilli)),2*ones(sum(receptor1.microvilli)),1) == [0,2,0,0,0,1,1,0]
@test filterphotons!(zeros(Int,8),receptor1,ones(sum(receptor1.microvilli)),2*ones(sum(receptor1.microvilli)),1) == [0,2,0,0,0,1,1,0]
@test filterphotons(receptor1,zeros(sum(receptor1.microvilli)),zeros(sum(receptor1.microvilli)),0) == vec(sum(receptor1.microvilli,2))

licvals1 = zeros(8) ; licvals1[1:length(bmp1)] = 2*bmp1 ; licvals1[7:8] += bmp1[1:2] ; licvals1[8] += 2*bmp1[1]
@test_approx_eq lightinducedcurrent(receptor1,zeros(sum(receptor1.microvilli)),zeros(sum(receptor1.microvilli)),bmp1) licvals1
@test_approx_eq lightinducedcurrent!(zeros(8),receptor1,zeros(sum(receptor1.microvilli)),zeros(sum(receptor1.microvilli)),bmp1) licvals1

#####################################################################################################################################

for args in [((:latency,:deterministic),GMHExamples.LatencyCalculation,:deterministic,(:latencyvalue,),(:uniform,)),
             ((:latency,:lognormal),GMHExamples.LatencyCalculation,:lognormal,(:latencylocation,:latencyscale),(:uniform,:uniform)),
             ((:refractory,:deterministic),GMHExamples.RefractoryCalculation,:deterministic,(:refractoryvalue,),(:uniform,)),
             ((:refractory,:lognormal),GMHExamples.RefractoryCalculation,:lognormal,(:refractorylocation,:refractoryscale),(:uniform,:uniform)),
             ((:bump,:fixed),GMHExamples.BumpShape,:fixed,(),()),
             ((:bump,:sample),GMHExamples.BumpShape,:sample,(:bumpamplitude,:bumpshape,:bumpscale),(:uniform,:lognormal,:lognormal))]
    t = trait(args[1]...)
    @test typeof(t) == args[2]
    @test traitvalue(t) == args[3]
    @test traittype(t) == Val{args[3]}
    @test GMHExamples.paramkeys(t) == args[4]
    @test GMHExamples.parampriors(t) == args[5]
end

b1 = fixedbumpshape()
b2 = fixedbumpshape(amplitude =3.0f0)
b3 = fixedbumpshape(shape =2.0f0)
b4 = fixedbumpshape(scale =1.5f0)
b5 = fixedbumpshape(cutoff =0.05f0)

for args in [(b1,4.0f0,3.0f0,2.5f0,0.025f0),
             (b2,3.0f0,3.0f0,2.5f0,0.025f0),
             (b3,4.0f0,2.0f0,2.5f0,0.025f0),
             (b4,4.0f0,3.0f0,1.5f0,0.025f0),
             (b5,4.0f0,3.0f0,2.5f0,0.05f0)]
    @test args[1].amplitude == args[2]
    @test args[1].shape == args[3]
    @test args[1].scale == args[4]
    @test args[1].cutoff == args[5]
end

p1 = policy(:photoreceptor)
p2 = policy(:photoreceptor,latency =:deterministic)
p3 = policy(:photoreceptor,refractory =:deterministic)
p4 = policy(:photoreceptor,bump =:sample)

for args in [(p1,:lognormal,:lognormal,:fixed,fixedbumpshape()),
             (p2,:deterministic,:lognormal,:fixed,fixedbumpshape()),
             (p3,:lognormal,:deterministic,:fixed,fixedbumpshape()),
             (p4,:lognormal,:lognormal,:sample,fixedbumpshape())]
    @test traitvalue(args[1].latency) == args[2]
    @test traitvalue(args[1].refractory) == args[3]
    @test traitvalue(args[1].bump) == args[4]
    @test args[1].fixedbump.amplitude == args[5].amplitude
    @test args[1].fixedbump.shape == args[5].shape
    @test args[1].fixedbump.scale == args[5].scale
    @test args[1].fixedbump.cutoff == args[5].cutoff
end

params1 = parameters(:photoreceptor,p1,(2.0,4.0),(0.2,0.8),(4.0,6.0),(1.5,2.5))
params2 = parameters(:photoreceptor,p2,(10.0,30.0),(4.0,6.0),(1.5,2.5))
params3 = parameters(:photoreceptor,p3,(2.0,4.0),(0.2,0.8),(50.0,150.0))
params4 = parameters(:photoreceptor,p4,(2.0,4.0),(0.2,0.8),(4.0,6.0),(1.5,2.5),(3.0,5.0),(log(3.0),0.3),(log(2.5),0.3))

for args in [(params1,[ParameterUnivariate,ParameterUnivariate,ParameterUnivariate,ParameterUnivariate]),
             (params2,[ParameterUnivariate,ParameterDefault,ParameterUnivariate,ParameterUnivariate]),
             (params3,[ParameterUnivariate,ParameterUnivariate,ParameterUnivariate,ParameterDefault]),
             (params4,repmat([ParameterUnivariate],7,1))]
    length(args[1]) == length(args[2])
    for i=1:length(args[1])
        typeof(args[1][i]) <: args[2][i]
    end
end

current = data(:array,1:5,rand(Float32,5))

photons = photonsequence(rphotons,100)
current = data(:array,0.001f0(1:100),zeros(Float32,100))

m1 = model(:photoreceptor,params1,photons,current,[1.0],10,p1)
m2 = model(:photoreceptor,params4,photons,current,[1.0],10,p4)

println("=======================")
show(m1)
show(m2)
println("=======================")

GMHExamples._latency!(Val{:deterministic},m1,[10.0,0.0,0.0,0.0])
GMHExamples._refractory!(Val{:deterministic},m1,[0.0,0.0,100.0,0.0])

@test m1.latency == 10.0f0*ones(Float32,length(m1.latency))
@test m1.refractory == 100.0f0*ones(Float32,length(m1.refractory))

latvals = zeros(Float32,length(m1.latency))
refvals = zeros(Float32,length(m1.refractory))
srand(467) ; GMHExamples._latency!(Val{:lognormal},m1,[2.5,0.5,0.0,0.0]) ; GMHExamples._refractory!(Val{:lognormal},m1,[0.0,0.0,5.0,0.7])
srand(467) ; GMHExamples.randlognormal!(latvals,2.5f0,0.5f0) ; GMHExamples.randlognormal!(refvals,5.0f0,0.7f0)
@test m1.latency == latvals
@test m1.refractory == refvals

@test GMHExamples._bumpvals(traittype(m1.policy.bump),m1,zeros(4)) == m1.fixedbumpvals
@test isempty(m2.fixedbumpvals) && GMHExamples._bumpvals(traittype(m2.policy.bump),m2,ones(7)) == bump(1.0f0,1.0f0,1.0f0,0.025f0)

GMHExamples._setval!(m1.filteredphotons,1)
GMHExamples._setval!(m1.refractimes,2)
GMHExamples._setval!(m1.lightinducedcurrent,3.0f0)

GMHExamples._resettemp(m1)

@test m1.filteredphotons == zeros(m1.filteredphotons)
@test m1.refractimes == zeros(m1.refractimes)
@test m1.lightinducedcurrent == zeros(m1.lightinducedcurrent)

evaluate!(m1,zeros(4))
evaluate!(m2,zeros(7))

@test size(evaluate!(m1,zeros(4),10)) == (100,10)
@test size(evaluate!(m2,zeros(7),10)) == (100,10)

println()
println("  * photoreceptor.jl tests completed *")



