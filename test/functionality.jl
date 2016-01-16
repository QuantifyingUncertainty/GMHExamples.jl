include("util.jl")

functionalitytests = [
  "plot_springmass",
  "plot_fitzhughnagumo",
  "plot_photoreceptor"
  ]

println("============================")
println("Running functionality tests:")
println("============================")

for t in functionalitytests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("functionality/",tfile))
end

