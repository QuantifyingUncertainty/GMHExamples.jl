#########################################
###Remove duplicates from a sorted vector
#########################################
function uniquesorted(in::AbstractVector)
  out = similar(in)
  c::Int = 0
  if length(in) > 0
    out[c+=1] = in[1]
    @inbounds for i=2:length(in)
      if in[i]!=in[i-1]
        out[c+=1]=in[i]
      end
    end
  end
  out[1:c]
end

############################################################
###Calculate a log-normal value for given location and scale
############################################################
@inline function randlognormal!{T<:AbstractFloat}(v::Vector{T},μ::T,σ::T)
    #tested different variations and the way this is implemented runs fastest
    #there is about 30% benefit in speed of using Float32 over Float64 for large vectors
    l = length(v)
    @simd for i=1:l
        @inbounds v[i] = T(randn())
    end
    @simd for i=1:l
        @inbounds v[i] = μ + σ*v[i]
    end
    @simd for i=1:l
        @inbounds v[i] = exp(v[i])
    end
    v
end

