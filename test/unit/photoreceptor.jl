nsteps = 5
nvilli = 8
rphotons = 0:10

###Test creation of photonsequence
p = photonsequence(rphotons,nsteps)
@test length(p) == nsteps
@test all(minimum(rphotons) .<= p .<= maximum(rphotons))

###Test utility functions
@test GMHExamples.uniquesorted(Int[]) == Int[]
@test GMHExamples.uniquesorted([2]) == [2]
@test GMHExamples.uniquesorted([1,1,1]) == [1]
@test GMHExamples.uniquesorted([1,2,2,3,3,3,8]) == [1,2,3,8]
@test GMHExamples.uniquesorted([1,4,7]) == [1,4,7]

###Test microvilli implementations
impl = [:full,:row,:column,:cell]
@test_throws ErrorException distributephotons(nvilli,p;impl=:dummy)
receptors = Dict{Symbol,Any}()
for i in impl
  srand(0)
  receptors[i] = distributephotons(nvilli,p;impl=i)
  @test numvilli(receptors[i]) == nvilli
  @test numsteps(receptors[i]) == nsteps
end

@test typeof(receptors[:full]) == GMHExamples.PhotoReceptorDense
@test typeof(receptors[:row]) == GMHExamples.PhotoReceptorRow
@test typeof(receptors[:column]) == GMHExamples.PhotoReceptorColumn
@test typeof(receptors[:cell]) == GMHExamples.PhotoReceptorCell{UInt16}
@test typeof(distributephotons(nvilli,p)) == GMHExamples.PhotoReceptorCell{UInt16}

for t=1:nsteps
  @test GMHExamples.microvilliwithphotons(receptors[:full],t) == GMHExamples.microvilliwithphotons(receptors[:row],t)
  @test GMHExamples.microvilliwithphotons(receptors[:full],t) == unique(GMHExamples.microvilliwithphotons(receptors[:cell],t))
end

for i=1:nvilli
  @test GMHExamples.photonimpacts(receptors[:full],i) == GMHExamples.photonimpacts(receptors[:column],i)
end

@test_throws ErrorException GMHExamples.microvilliwithphotons(receptors[:column],1)
@test_throws ErrorException GMHExamples.photonimpacts(receptors[:row],1)
@test_throws ErrorException GMHExamples.photonimpacts(receptors[:cell],1)

###Test current integration functions
@test GMHExamples._bumpvalmaximum(0.0,0.0,0.0) == 0.0
@test GMHExamples._bumpvalmaximum(2.0,3.0,2.5) == 7.5
@test GMHExamples._bumpvalmaximum(2.0f0,3.0f0,2.5f0) == 7.5f0 #Float32
@test GMHExamples._bumpval(1.0,1.0,1.0,0.0) == 0.0
@test GMHExamples._bumpval(2.0,3.0,1.5,0.0) == 0.0
@test GMHExamples._bumpval(2.0,3.0,1.5,4.5) == 2.0
@test GMHExamples._bumpval(4.0f0,3.0f0,2.5f0,7.5f0) == 4.0f0 #Float32
@test bump!(4.0f0,3.0f0,2.5f0,Array{Float32}(30)) == bump(4.0f0,3.0f0,2.5f0,30)
@test bump!(4.0,3.0,2.5,Array{Float64}(30)) == bump(4.0,3.0,2.5)
@test length(bump(0.0,0.0,0.0)) == 0

bmp1 = bump(1.0,1.0,1.0)
@test maximum(bmp1) == 1.0 && GMHExamples._bumpvalmaximum(1.0,1.0,1.0) == 1.0 && length(bmp1) == 6

ph = zeros(Int,8,7) ; mc = zeros(8,7)
ph[1,1] = 1 ; mc[1:6,1] = bmp1
ph[2,2] = 1 ; mc[2:7,2] = bmp1
ph[4,3] = 1 ; mc[4:8,3] = bmp1[1:5]
ph[1,4] = 2 ; mc[:,4] = 2*mc[:,1]
ph[1:2,5] = [1,1] ; mc[:,5] = mc[:,1] + mc[:,2]
ph[[2,4],6] = [1,1] ; mc[:,6] = mc[:,2] + mc[:,3]
ph[[1,2,4],7] = [2,1,1] ; mc[:,7] = mc[:,4] + mc[:,6]

for j=1:size(ph,2)
  @test_approx_eq mc[:,j] macrocurrent!(ph[:,j],bmp1,zeros(8))
  @test_approx_eq mc[:,j] macrocurrent(ph[:,j],bmp1)
end
@test_throws AssertionError macrocurrent!(ph[:,1],bmp1,zeros(10))

println()
println("  * photoreceptor.jl tests completed *")



