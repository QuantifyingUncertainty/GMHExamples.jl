##############
println("Loading photon sequence from file")
gc() ; photons = photonsequence("photo/data/naturallight.jld")
for i=1:3
    gc() ; @time photons = photonsequence("photo/data/naturallight.jld")
end
println("Loading light-induced current from file")
gc() ; lic = lightinducedcurrent("photo/data/naturallight.jld")
for i=1:3
    gc() ; @time lic = lightinducedcurrent("photo/data/naturallight.jld")
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
nsteps = numvalues(photons)
receptorimpl = [:full,:row,:cell]
receptors = Dict{Symbol,Any}()

for args in [(:full,(UInt8,)),(:row,(UInt16,Int)),(:cell,(UInt16,))]
  receptors[args[1]] = photoreceptor(args[1],datavalues(photons),nvilli,args[2]...)
  println("Creating AbstractPhotoReceptor of type ",args[1])
  for j=1:3
    receptors[args[1]] = []
    gc()
    @time receptors[args[1]] = photoreceptor(args[1],datavalues(photons),nvilli,args[2]...)
  end
end

@inline function testmicrovilliaccess(receptor::GMHExamples.AbstractPhotoReceptor)
    a = 0
    for t=1:numtimesteps(receptor)
        microvilli = GMHExamples.microvilliwithphotons(receptor,t)
        for i=1:length(microvilli)
            a += microvilli[i]
        end
    end
    a
end

for i in [:full,:row,:cell]
  testmicrovilliaccess(receptors[i])
  println("Testing microvilliwithphotons access of type ",i)
  for j=1:3
    gc() ; @time testmicrovilliaccess(receptors[i])
  end
end

gc()
for i in [:full,:row,:cell]
  println("Memory taken by Microvilli object of type ",i,": ",objectsizetostr(receptors[i]))
end

gc()
for T in [Float32,Float64]
    println("=================")
    println("Testing ",T)
    println("=================")
    ###Testing of bump functions
    timepoints = collect(T(1.0):T(1.0):T(30.0))
    shape = T[4.0,3.0,2.5]
    cutoff = T(0.025)
    bumpvals = zeros(T,30)
    bump!(bumpvals,timepoints,shape[1],shape[2],shape[3])
    bump(timepoints,shape[1],shape[2],shape[3])
    bump(shape[1],shape[2],shape[3],cutoff)
    println("Testing bump!{",T,"}")
    for j=1:3
        @time bumpvals = bump!(bumpvals,timepoints,shape[1],shape[2],shape[3])
    end
    println("Testing bump{",T,"} with given length")
    for j=1:3
        @time bumpvals = bump(timepoints,shape[1],shape[2],shape[3])
    end
    println("Testing bump{",T,"} with tolerance")
    for j=1:3
        @time bumpvals = bump(shape[1],shape[2],shape[3],cutoff)
    end
    bumplength = length(bumpvals)
    ###Testing of macrocurrent functions
    current = zeros(T,nsteps)
    macrocurrent!(current,datavalues(photons),bumpvals)
    macrocurrent(datavalues(photons),bumpvals)
    gc()
    println("Testing macrocurrent!{",T,"}")
    for j=1:3
        current = zeros(T,nsteps)
        @time macrocurrent!(current,datavalues(photons),bumpvals)
    end
    gc()
    println("Testing macrocurrent{",T,"}")
    for j=1:3
        @time current = macrocurrent(datavalues(photons),bumpvals)
    end

    latparas = T[2.7137,0.4568]
    refparas = T[log(150),log(1.7)]
    latency = zeros(T,sum(datavalues(photons)))
    refractory = zeros(T,sum(datavalues(photons)))
    GMHExamples.randlognormal!(latency,latparas[1],latparas[2])
    GMHExamples.randlognormal!(refractory,refparas[1],refparas[2])
    println("Generating latency and refractory values")
    for j=1:3
        gc() ; @time (GMHExamples.randlognormal!(latency,latparas[1],latparas[2]) ; GMHExamples.randlognormal!(refractory,refparas[1],refparas[2]))
    end

    for impl in [:cell]
        GMHExamples._filterphotons!(zeros(Int,nsteps),zeros(Int,nvilli),receptors[impl],latency,refractory,bumplength)
        println("Testing _filterphotons!{$(T)} for $(impl)")
        for i=1:5
            gc()
            fphotons = zeros(Int,nsteps)
            rtimes = zeros(Int,nvilli)
            @time GMHExamples._filterphotons!(fphotons,rtimes,receptors[impl],latency,refractory,bumplength)
        end
        filterphotons!(zeros(Int,nsteps),receptors[impl],latency,refractory,bumplength)
        println("Testing filterphotons!{$(T)} for $(impl)")
        for i=1:5
            gc()
            fphotons = zeros(Int,nsteps)
            @time filterphotons!(fphotons,receptors[impl],latency,refractory,bumplength)
        end
        filterphotons(receptors[impl],latency,refractory,bumplength)
        println("Testing filterphotons{$(T)} for $(impl)")
        for i=1:5
            gc() ; @time filterphotons(receptors[impl],latency,refractory,bumplength)
        end
        ######################################################################################
        GMHExamples._lightinducedcurrent!(zeros(T,nsteps),zeros(Int,nsteps),zeros(Int,nvilli),receptors[impl],latency,refractory,bumpvals)
        println("Testing _lightinducedcurrent!{$(T)} for $(impl)")
        for i=1:5
            gc()
            current = zeros(T,nsteps)
            fphotons = zeros(Int,nsteps)
            rtimes = zeros(Int,nvilli)
            @time GMHExamples._lightinducedcurrent!(current,fphotons,rtimes,receptors[impl],latency,refractory,bumpvals)
        end
        lightinducedcurrent!(zeros(T,nsteps),receptors[impl],latency,refractory,bumpvals)
        println("Testing lightinducedcurrent!{$(T)} for $(impl)")
        for i=1:5
            gc()
            current = zeros(T,nsteps)
            @time lightinducedcurrent!(current,receptors[impl],latency,refractory,bumpvals)
        end
        lightinducedcurrent(receptors[impl],latency,refractory,bumpvals)
        println("Testing lightinducedcurrrent{$(T)} for $(impl)")
        for i=1:5
            gc() ; @time lightinducedcurrent(receptors[impl],latency,refractory,bumpvals)
        end
    end

    licit = 10
    licvec = Array(T,numvalues(photons),licit)

    println("$licit iterations of lightinducedcurrent for :cell")
    for i=1:licit
        gc() ; @time (GMHExamples.randlognormal!(latency,latparas[1],latparas[2]) ; GMHExamples.randlognormal!(refractory,refparas[1],refparas[2]) ; licvec[:,i] = lightinducedcurrent(receptors[:cell],latency,refractory,bumpvals))
    end

    println("100 iterations of lightinducedcurrent! for :cell")
    gc()
    @time for i=1:100
        GMHExamples.randlognormal!(latency,latparas[1],latparas[2])
        GMHExamples.randlognormal!(refractory,refparas[1],refparas[2])
        lightinducedcurrent!(current,receptors[:cell],latency,refractory,bumpvals)
    end
end

current = data(:array,0.001f0(1:nsteps),zeros(Float32,nsteps))
p1 = policy(:photoreceptor)
p2 = policy(:photoreceptor,bump =:sample)

params1 = parameters(:photoreceptor,p1,(2.0,4.0),(0.2,0.8),(4.0,6.0),(1.5,2.5))
params2 = parameters(:photoreceptor,p2,(2.0,4.0),(0.2,0.8),(4.0,6.0),(1.5,2.5),(3.0,5.0),(log(3.0),0.3),(log(2.5),0.3))

m1 = model(:photoreceptor,params1,photons,current,[1.0],30000,p1)
m2 = model(:photoreceptor,params2,photons,current,[1.0],30000,p2)

println("Evaluating the model with latency and refractory parameters")
evaluate!(m1,[2.7,0.5,5.0,2.0])
for i=1:5
    gc() ; @time evaluate!(m1,[2.7,0.5,5.0,2.0])
end

println("Evaluating the model with latency, refractory and bump parameters")
evaluate!(m2,[2.7,0.5,5.0,2.0,4.0,3.0,2.5])
for i=1:5
    gc() ; @time evaluate!(m2,[2.7,0.5,5.0,2.0,4.0,3.0,2.5])
end

println("Evaluating the model with latency and refractory parameters 1000 times")
gc()
@time for i=1:1000
    evaluate!(m1,[2.7,0.5,5.0,2.0,4.0,3.0,2.5])
end

results1 = evaluate!(m1,[2.7,0.5,5.0,2.0],1)
println("Evaluating the model 1000 times with latency and refractory parameters, storing the results")
@time results1 = evaluate!(m1,[2.7,0.5,5.0,2.0],1000)

println("Evaluating the model 1000 times with latency, refractory and bump parameters, storing the results")
@time results2 = evaluate!(m2,[2.7,0.5,5.0,2.0,4.0,3.0,2.5],1000)
