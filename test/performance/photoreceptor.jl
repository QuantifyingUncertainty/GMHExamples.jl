##############
println("Loading photon sequence from file")
gc() ; p = photonsequence("photo/data/naturallight.jld")
for i=1:3
  gc() ; @time p = photonsequence("photo/data/naturallight.jld")
end
###############
testunique() = (for j=1:2000 unique(sort!(rand!(Array{Int}(150),1:30000))) end ; nothing)
testuniquesorted() = (for i=1:2000 GMHExamples.uniquesorted(sort!(rand!(Array{Int}(150),1:30000))) end ; nothing)
testunique()
testuniquesorted()
gc()
println("Testing unique")
for i=1:3
  gc()
  @time testunique()
end
println("Testing uniquesorted")
for i=1:3
  gc()
  @time testuniquesorted()
end
###############
gc()
nvilli = 30000
nsteps = length(p)
receptorimpl = [:full,:row,:column,:cell]
receptors = Dict{Symbol,Any}()

for i in receptorimpl
  receptors[i] = distributephotons(nvilli,p;impl=i)
  println("Creating PhotoReceptor of type ",i)
  for j=1:3
    receptors[i] = []
    gc()
    @time receptors[i] = distributephotons(nvilli,p;impl=i)
  end
end

testmicrovilliaccess(receptor::PhotoReceptor) = (a = 0 ; for t=1:numsteps(receptor) microvilli = GMHExamples.microvilliwithphotons(receptor,t) ; for i=1:length(microvilli) a += microvilli[i] end end ; a)
testimpactaccess(receptor::PhotoReceptor) = (a = 0 ; for i=1:numvilli(receptor) impacts = GMHExamples.photonimpacts(receptor,i) ; for t=1:length(impacts) a += impacts[t] end end ; a)

for i in [:full,:row,:cell]
  testmicrovilliaccess(receptors[i])
  println("Testing microvilliwithphotons access of type ",i)
  for j=1:3
    gc() ; @time testmicrovilliaccess(receptors[i])
  end
end

for i in [:full,:column]
  testimpactaccess(receptors[i])
  println("Testing photonimpacts access of type ",i)
  for j=1:3
    gc() ; @time testimpactaccess(receptors[i])
  end
end

gc()
for i in receptorimpl
  println("Memory taken by Microvilli object of type ",i,": ",objectsizetostr(receptors[i]))
end

gc()
for T in [Float32,Float64]
  println("=================")
  println("Testing ",T)
  println("=================")
  ###Testing of bump functions
  shape = T[4.0,3.0,2.5]
  bumpvals = zeros(T,30)
  bump!(shape[1],shape[2],shape[3],bumpvals)
  bump(shape[1],shape[2],shape[3],30)
  bump(shape[1],shape[2],shape[3])
  println("Testing bump!{",T,"}")
  for j=1:3
    @time bumpvals = bump!(shape[1],shape[2],shape[3],bumpvals)
  end
  println("Testing bump{",T,"} with given length")
  for j=1:3
    @time bumpvals = bump(shape[1],shape[2],shape[3],30)
  end
  println("Testing bump{",T,"} with tolerance")
  for j=1:3
    @time bumpvals = bump(shape[1],shape[2],shape[3])
  end
  ###Testing of macrocurrent functions
  current = zeros(T,nsteps)
  macrocurrent!(p,bumpvals,current)
  macrocurrent(p,bumpvals)
  gc()
  println("Testing macrocurrent!{",T,"}")
  for j=1:3
    current = zeros(T,nsteps)
    @time macrocurrent!(p,bumpvals,current)
  end
  gc()
  println("Testing macrocurrent{",T,"}")
  for j=1:3
    @time current = macrocurrent(p,bumpvals)
  end
end

# latparas = [2.7137,0.4568]
# refparas = [log(150),log(1.7)]
# latpdf = GMHExamples.latencypdf(latparas[1],latparas[2])
# refpdf = GMHExamples.refractorypdf(refparas[1],refparas[2])

# rt1 = lic(receptors[:cell],bumpvals,10,100)
# rt2 = lic(receptors[:cell],bumpvals,latpdf,refpdf)
# rt3 = lic(receptors[:row],bumpvals,latpdf,refpdf)
# println("Testing deterministic lic")
# for i=1:5
#   gc() ; @time rt1 = lic(receptors[:cell],bumpvals,10,100)
# end
# println("Testing stochastic lic with :cell")
# for i=1:5
#   gc() ; @time rt2 = lic(receptors[:cell],bumpvals,latpdf,refpdf)
# end
# println("Testing stochastic lic with :row")
# for i=1:5
#   gc() ; @time rt3 = lic(receptors[:row],bumpvals,latpdf,refpdf)
# end

# rt1vec = Array(Float64,length(p),10)
# rt2vec = Array(Float64,length(p),10)

# l1vec[:,1] = lightinducedcurrent(d1,shapeparas,10,100)
# l2vec[:,1] = lightinducedcurrent(d1,shapeparas,latencyparas,refracparas)
# println("Calculating light-induced current (deterministic latency and refractory parameters)")
# for i=1:10
#   gc() ; @time l1vec[:,i] = lightinducedcurrent(d1,shapeparas,10,100)
# end

# println("Calculating light-induced current (stochastic latency and refractory parameters)")
# for i=1:10
#   gc() ; @time l2vec[:,i] = lightinducedcurrent(d1,shapeparas,latencyparas,refracparas)
# end

# println("100 Iterations of lightinducedcurrent (deterministic)")
# l3vec = Array(Float64,length(p),100)
# @time for i=1:10
#   l3vec[:,i] = lightinducedcurrent(d1,shapeparas,10,100)
# end

# println("100 Iterations of lightinducedcurrent (stochastic)")
# l4vec = Array(Float64,length(p),100)
# @time for i=1:10
#   l4vec[:,i] = lightinducedcurrent(d1,shapeparas,latencyparas,refracparas)
# end





