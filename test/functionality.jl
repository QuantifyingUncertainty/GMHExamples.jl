functionalitytests = [
  "plot_springmass",
  "plot_fitzhughnagumo"
  ]

println("============================")
println("Running functionality tests:")
println("============================")

for t in functionalitytests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("functionality/",tfile))
end

