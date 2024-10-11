using TOML

svp = include("/home/brian/Code/SVP/src/SVP.jl")

inp = TOML.parsefile(ARGS[1])

for k in keys(inp)
  println(k)
  re = svp.readHess(inp[k]["reactant"])
  ts = svp.readHess(inp[k]["ts"])
  pr = svp.readHess(inp[k]["product"])
  
  if haskey(inp[k], "molecule") && haskey(inp[k], "surface")
    mol  = inp[k]["molecule"]
    surf = inp[k]["surface"]
    svp.printSVP(re, ts, pr, mol=mol, surf=surf) |> println
  else
    svp.printSVP(re, ts, pr) |> println
  end

  println("End $(k)\n")
end
