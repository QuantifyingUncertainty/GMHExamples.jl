################
###Photoreceptor
################
abstract PhotoReceptor

numsteps(r::PhotoReceptor) = r.numsteps
numvilli(r::PhotoReceptor) = r.numvilli
microvilliwithphotons(r::PhotoReceptor,t::Int) = error("Not implemented because microvilliwithphotons is not efficient for PhotoReceptor of type ",typeof(r))
photonimpacts(r::PhotoReceptor,i::Int) = error("Not implemented because photonimpacts is not efficient for PhotoReceptor of type ",typeof(r))

#################################
###Photons stored in dense matrix
#################################

type PhotoReceptorDense <: PhotoReceptor
  numsteps::Int
  numvilli::Int
  microvilli::AbstractMatrix{UInt8}
end

function PhotoReceptorDense(nvilli::Int,photons::AbstractVector{Int})
  nsteps = length(photons)
  m = zeros(UInt8,nsteps,nvilli)
  @inbounds for i = 1:nsteps
    @inbounds for j in rand(1:nvilli,photons[i])
      m[i,j] += 0x01
    end
  end
  PhotoReceptorDense(nsteps,nvilli,m)
end

microvilliwithphotons(r::PhotoReceptorDense,t::Int) = find(r.microvilli[t,:])
photonimpacts(r::PhotoReceptorDense,i::Int) = find(r.microvilli[:,i])

###########################################################
###Photons stored in sparse matrix, one microvillus per row
###########################################################

type PhotoReceptorRow <: PhotoReceptor
  numsteps::Int
  numvilli::Int
  microvilli::AbstractMatrix{UInt8}
end

function PhotoReceptorRow(nvilli::Int,photons::AbstractVector{Int})
  nsteps = length(photons)
  e = cumsum(photons) ; b = e - photons + 1 ; nphotons = e[end]
  tp = Array{Int}(nphotons)
  @inbounds for i=1:nsteps
    tp[b[i]:e[i]] = i
  end
  PhotoReceptorRow(nsteps,nvilli,sparse(rand(1:nvilli,nphotons),tp,ones(UInt8,nphotons),nvilli,nsteps))
end

microvilliwithphotons(r::PhotoReceptorRow,t::Int) = sub(rowvals(r.microvilli),nzrange(r.microvilli,t))

##############################################################
###Photons stored in sparse matrix, one microvillus per column
##############################################################

type PhotoReceptorColumn <: PhotoReceptor
  numsteps::Int
  numvilli::Int
  microvilli::AbstractMatrix{UInt8}
end

function PhotoReceptorColumn(nvilli::Int,photons::AbstractVector{Int})
  nsteps = length(photons)
  e = cumsum(photons) ; b = e - photons + 1 ; nphotons = e[end]
  tp = Array{Int}(nphotons)
  @inbounds for i=1:nsteps
    tp[b[i]:e[i]] = i
  end
  PhotoReceptorColumn(nsteps,nvilli,sparse(tp,rand(1:nvilli,nphotons),ones(UInt8,nphotons),nsteps,nvilli))
end

photonimpacts(r::PhotoReceptorColumn,i::Int) = sub(rowvals(r.microvilli),nzrange(r.microvilli,i))

################################
###Photons stored in cell matrix
################################

type PhotoReceptorCell{T<:Integer} <: PhotoReceptor
  numsteps::Int
  numvilli::Int
  microvilli::Vector{Vector{T}}
end

PhotoReceptorCell(nvilli::Int,photons::AbstractVector{Int}) = PhotoReceptorCell(length(photons),nvilli,[sort!(rand!(Vector{UInt16}(p),1:nvilli)) for p in photons])

microvilliwithphotons(r::PhotoReceptorCell,t::Int) = r.microvilli[t]

##################################
###Factory function for microvilli
##################################

function distributephotons(nvilli::Int,photons::AbstractVector{Int};impl::Symbol = :cell)
  if impl == :full
    PhotoReceptorDense(nvilli,photons)
  elseif impl == :row
    PhotoReceptorRow(nvilli,photons)
  elseif impl == :column
    PhotoReceptorColumn(nvilli,photons)
  elseif impl == :cell
    PhotoReceptorCell(nvilli,photons)
  else
    error("Unknown PhotoReceptor implementation $impl")
  end
end
