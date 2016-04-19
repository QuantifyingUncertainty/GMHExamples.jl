srand(0)

unittests = [
  "photoreceptor"]

println("===================")
println("Running unit tests:")
println("===================")

for t in unittests
  tfile = t*".jl"
  println("  * $(tfile) *")
  include(string("unit/",tfile))
  println()
  println()
end

