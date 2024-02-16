using TOML

svp = include("/home/brian/Code/SVP/src/SVP.jl")

inp = TOML.parsefile(ARGS[1])

for k in keys(inp)
  println(k)
  re = svp.readHess(inp[k]["reactant"])
  ts = svp.readHess(inp[k]["ts"])
  pr = svp.readHess(inp[k]["product"])
  
  svp.printSVP(re, ts, pr) |> println

  println("End $(k)\n")
end
