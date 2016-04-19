################
###Photoreceptor
################
abstract AbstractPhotoReceptor

numtimesteps(r::AbstractPhotoReceptor) = r.numsteps
numvilli(r::AbstractPhotoReceptor) = r.numvilli
numphotons(r::AbstractPhotoReceptor) = r.numphotons

photoreceptor(s::Symbol,photons::Vector,nvilli::Integer,args...) = _photoreceptor(Val{s},photons,nvilli,args...)

###############################################
###Photons stored in a matrix (dense or sparse)
###############################################

type PhotoReceptorMatrix{T<:Integer,A<:AbstractMatrix} <: AbstractPhotoReceptor
    microvilli::A
    numsteps::Int
    numvilli::Int
    numphotons::Int
    PhotoReceptorMatrix(mv::AbstractMatrix{T},ns::Int,nv::Int,np::Int) = new(mv,ns,nv,np)
end

### Dense Matrix factory functions

function _photoreceptor{T<:Integer}(::Type{Val{:full}},photons::Vector,nvilli::Integer,::Type{T})
    nsteps = length(photons)
    nphotons = sum(photons)
    m = zeros(T,nsteps,nvilli)
    for i = 1:nsteps
        r = rand(1:nvilli,photons[i])
        for j in r
            @inbounds m[i,j] += one(T)
        end
    end
    PhotoReceptorMatrix{T,typeof(m)}(m,nsteps,nvilli,nphotons)
end

### Sparse Matrix factory functions (row version)

@inline function _photontimeindex{I<:Integer}(photons::Vector,::Type{I})
    e = cumsum(photons)
    photontime = Vector{I}(e[end])
    for i=one(I):I(length(photons))
        @simd for j=e[i]-photons[i]+1:e[i]
            @inbounds photontime[j] = i
        end
    end
    photontime
end

@inline _randvilliindex{I<:Integer}(nvilli::Int,nphotons::Int,::Type{I}) = rand(one(I):I(nvilli),nphotons)

function _photoreceptor{T<:Integer,I<:Integer}(::Type{Val{:row}},photons::Vector,nvilli::Integer,::Type{T},::Type{I})
    pti = _photontimeindex(photons,I)
    nphotons = length(pti)
    nsteps = length(photons)
    m = sparse(_randvilliindex(nvilli,nphotons,I),pti,ones(T,nphotons),nvilli,nsteps)
    PhotoReceptorMatrix{T,typeof(m)}(m,nsteps,nvilli,nphotons)
end

@inline microvilliwithphotons{I<:Integer,M<:DenseMatrix}(r::PhotoReceptorMatrix{I,M},timepoint::Integer) = find(r.microvilli[timepoint,:])
@inline microvilliwithphotons{I<:Integer,S<:SparseMatrixCSC}(r::PhotoReceptorMatrix{I,S},timepoint::Integer) = sub(rowvals(r.microvilli),nzrange(r.microvilli,timepoint))

################################
###Photons stored in cell matrix
################################

type PhotoReceptorCell{T<:Integer} <: AbstractPhotoReceptor
    microvilli::Vector{Vector{T}}
    numsteps::Int
    numvilli::Int
    numphotons::Int
    PhotoReceptorCell(mv::Vector{Vector{T}},ns::Int,nv::Int,np::Int) = new(mv,ns,nv,np)
end

function _photoreceptor{T<:Integer}(::Type{Val{:cell}},photons::Vector,nvilli::Integer,::Type{T})
    PhotoReceptorCell{T}([sort!(rand!(Vector{T}(p),one(T):T(nvilli))) for p in photons],length(photons),nvilli,sum(photons))
end

@inline microvilliwithphotons(r::PhotoReceptorCell,timepoint::Integer) = r.microvilli[timepoint]

function show(io::IO,r::AbstractPhotoReceptor)
    println(io,"PhotoReceptor with $(r.numvilli) microvilli and $(r.numsteps) timesteps")
    println(io,"  Additional fields: :microvilli")
    nothing
end
