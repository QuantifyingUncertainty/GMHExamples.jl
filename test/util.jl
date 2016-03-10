using Base.Test
using GMHExamples
using GeneralizedMetropolisHastings

function objectsizetostr(a::Any)
  bytes = Base.summarysize(a)
  result = ""
  if bytes < 10*1024
    result = @sprintf("%6d bytes  ",bytes)
  elseif bytes < 1024*1024
    result = @sprintf("%5.3f kB",bytes/(1024))
  else
    result = @sprintf("%5.3f MB",bytes/1024/1024)
  end
  result
end

const LN2 = 0.6931471805599453
const EXPA = Int(2^52/LN2)
const EXPB = 1023*2^52
const EXPC = 60801
const EXPBC = EXPB - EXPC

expfast(x::Float64) = reinterpret(Float64,Int(EXPA*x) + EXPBC)
expfast(x::AbstractVector) = [expfast(i) for i in x]

expcf1(x::Float64,ei::Float64) = (1.0 + 2.0*x/(2.0 - x + 0.16666666666666666*x*x))*ei
expcf2(x::Float64,ei::Float64) = (1.0 + 2.0*x/(2.0 - x + 0.16666666666666666*x*x/(1.0 + 0.01666666666666666*x*x)))*ei
expcf3(x::Float64,ei::Float64) = (1.0 + 2.0*x/(2.0 - x + 0.16666666666666666*x*x/(1.0 + 0.01666666666666666*x*x/(1.0 + 0.007142857142857143*x*x))))*ei

expcf1(x::Float32,ei::Float32) = (1.0f0 + 2.0f0*x/(2.0f0 - x + 0.16666667f0*x*x))*ei
expcf2(x::Float32,ei::Float32) = (1.0f0 + 2.0f0*x/(2.0f0 - x + 0.16666667f0*x*x/(1.0f0 + 0.016666668f0*x*x)))*ei
expcf3(x::Float32,ei::Float32) = (1.0f0 + 2.0f0*x/(2.0f0 - x + 0.16666667f0*x*x/(1.0f0 + 0.016666668f0*x*x/(1.0f0 + 0.007142857f0*x*x))))*ei

testsum{T<:Real}(x::AbstractVector{T},μ::T,σ::T) = (a::T = zero(T) ; @simd for i=1:length(x)  @inbounds a+=(μ + σ*x[i]) end ; a)
testexp{T<:Real}(x::AbstractVector{T},μ::T,σ::T) = (a::T = zero(T) ; @simd for i=1:length(x)  @inbounds a+=exp(μ + σ*x[i]) end ; a)
testexpfast(x::AbstractVector{Float64},μ::Float64,σ::Float64) = (a = 0.0 ; @simd for i=1:length(x) @inbounds a+=expfast(μ + σ*x[i]) end ; a)

function testexpcf1{T<:Real}(x::AbstractVector{T},μ::T,σ::T)
  a::T = zero(T)
  y::T = zero(T)
  yr::T = zero(T)
  e = [exp(s) for s=T(0.0):T(10.0)]
  @simd for i=1:length(x)
    @inbounds y = (μ + σ*x[i])
    yr = y % 1.0
    @inbounds a += expcf1(yr,e[trunc(Int,y-yr)])
  end
  a
end

function testexpcf2{T<:Real}(x::AbstractVector{T},μ::T,σ::T)
  a::T = zero(T)
  y::T = zero(T)
  yr::T = zero(T)
  e = [exp(s) for s=T(0.0):T(10.0)]
  @simd for i=1:length(x)
    @inbounds y = (μ + σ*x[i])
    yr = y % 1.0
    @inbounds a += expcf2(yr,e[trunc(Int,y-yr)])
  end
  a
end
function testexpcf3{T<:Real}(x::AbstractVector{T},μ::T,σ::T)
  a::T = zero(T)
  y::T = zero(T)
  yr::T = zero(T)
  e = [exp(s) for s=T(0.0):T(10.0)]
  @simd for i=1:length(x)
    @inbounds y = (μ + σ*x[i])
    yr = y % 1.0
    @inbounds a += expcf3(yr,e[trunc(Int,y-yr)])
  end
  a
end

function testall{T<:Real}(x::AbstractVector{T},μ::T,σ::T)
  for f in [testsum,testexp,testexpcf1,testexpcf2,testexpcf3]
    f(x,μ,σ)
  end
  println("Test of type ",typeof(μ))
  gc()
  for f in [testsum,testexp,testexpcf1,testexpcf2,testexpcf3]
    gc()
    println(f)
    for i=1:5
      @time f(x,μ,σ)
    end
  end
end


@inline randnormal{T<:AbstractFloat}(μ::T,σ::T) = @fastmath μ + σ*T(randn())
@inline randlognormal{T<:AbstractFloat}(μ::T,σ::T) = @fastmath exp(μ + σ*T(randn()))

testrandnormal(μ::AbstractFloat,σ::AbstractFloat,l::Int) = @simd for i=1:l randnormal(μ,σ) end
testrandlognormal(μ::AbstractFloat,σ::AbstractFloat,l::Int) = @simd for i=1:l randlognormal(μ,σ) end

function testrandnormal{T<:AbstractFloat}(v::Vector{Int},μ::T,σ::T)
    @simd for i=1:length(v)
        @inbounds v[i] = trunc(Int,randnormal(μ,σ))
    end
    v
end

function testrandlognormal{T<:AbstractFloat}(v::Vector{Int},μ::T,σ::T)
    @simd for i=1:length(v)
        r = randlognormal(μ,σ)
        @inbounds v[i] = trunc(Int,r)
    end
    v
end

function testrandlognormal{T<:AbstractFloat}(v::Vector{T},μ::T,σ::T)
    @simd for i=1:length(v)
        @inbounds v[i] = T(randn())
    end
    @simd for i=1:length(v)
        @inbounds v[i] = μ + σ*v[i]
    end
    @simd for i=1:length(v)
        @inbounds v[i] = exp(v[i])
    end
end

function testrandlognormaltrunc{T<:AbstractFloat}(v::Vector{T},μ::T,σ::T)
    a = 0
    testrandlognormal(v,μ,σ)
    @simd for i=1:length(v)
        @inbounds a += trunc(Int,v[i])
    end
end

function testrandlognormal{T<:AbstractFloat}(iv::Vector{Int},v::Vector{T},μ::T,σ::T)
    @simd for i=1:length(v)
        @inbounds v[i] = T(randn())
    end
    @simd for i=1:length(v)
        @inbounds v[i] = μ + σ*v[i]
    end
    @simd for i=1:length(v)
        @inbounds v[i] = exp(v[i])
    end
    @simd for i=1:length(v)
        @inbounds iv[i] = trunc(Int,v[i])
    end
end
