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


randnormal(μ::Float64,σ::Float64) = μ + σ*randn()
randlognormal(μ::Float64,σ::Float64) = exp(μ + σ*randn())
randlognormalfast(μ::Float64,σ::Float64) = expfast(μ + σ*randn())

randlognormalcf0(μ::Float64,σ::Float64) = expcf3(μ + σ*randn())
randlognormalcf1(μ::Float64,σ::Float64) = (r::Float64 = μ + σ*randn() ; r2::Int = trunc(Int,r/LN2) ; (1<<r2)*expcf1(r-r2*LN2))
randlognormalcf2(μ::Float64,σ::Float64) = (r::Float64 = μ + σ*randn() ; r2::Int = trunc(Int,r/LN2) ; (1<<r2)*expcf2(r-r2*LN2))
randlognormalcf3(μ::Float64,σ::Float64) = (r::Float64 = μ + σ*randn() ; r2::Int = trunc(Int,r/LN2) ; (1<<r2)*expcf3(r-r2*LN2))

testrandnormal(μ::Float64,σ::Float64,l::Int) = for i=1:l randnormal(μ,σ) end
testrandlognormal(μ::Float64,σ::Float64,l::Int) = for i=1:l randlognormal(μ,σ) end
testrandlognormalfast(μ::Float64,σ::Float64,l::Int) = for i=1:l randlognormalfast(μ,σ) end
testrandlognormalcf0(μ::Float64,σ::Float64,l::Int) = for i=1:l randlognormalcf0(μ,σ) end
testrandlognormalcf1(μ::Float64,σ::Float64,l::Int) = for i=1:l randlognormalcf1(μ,σ) end
testrandlognormalcf2(μ::Float64,σ::Float64,l::Int) = for i=1:l randlognormalcf2(μ,σ) end
testrandlognormalcf3(μ::Float64,σ::Float64,l::Int) = for i=1:l randlognormalcf3(μ,σ) end


