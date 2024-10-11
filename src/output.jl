
function CoM(bdys)
  M = sum([i.m for i in bdys])
  r = sum([i.m*i.r for i in bdys])
  return r ./ M
end

function printSVP(reac, ts, prod; mol=(), surf=())
  io = IOBuffer()

  rc = let 
    i = findfirst(e -> e < 0, ts.freq)
    ts.modes[:, i]
  end

  n = size(reac.modes)[2]
  N = div(n, 3)
  c = sqrt(N)

  X = repeat([1.0, 0.0, 0.0], N) ./ c
  Y = repeat([0.0, 1.0, 0.0], N) ./ c
  Z = repeat([0.0, 0.0, 1.0], N) ./ c
  
  println(io, "Forward Direction")
  println(io, "-----------------")
  println(io, " Mode    SVP")
  
  for i in 1:n
    p = dot(reac.modes[:, i], rc) |> abs |> (x -> round(x, digits=5))
    println(io, "  $(i)      $(p)")
  end

  println(io, "\n")
  println(io, "Backword Direction")
  println(io, "------------------")
  println(io, " Mode    SVP")
  
  for i in 1:n
    p = dot(prod.modes[:, i], rc) |> abs |> (x -> round(x, digits=5))
    println(io, "  $(i)      $(p)")
  end

  println(io, "\n")
  println(io, "Translational")
  println(io, "-------------")
  println(io, " Axis    SVP")

  for a in [(X,"x"), (Y,"y"), (Z,"z")]
    p = dot(a[1], rc) |> abs |> (x -> round(x, digits=5))
    println(io, "  $(a[2])      $(p)")
  end

  if isempty(mol) || isempty(surf)
    return take!(io) |> String
  end

  println(io, "\n")
  println(io, "CoM Difference Component")
  println(io, "------------------------")

  molCoM  = CoM(prod.bdys[mol[1]:mol[2]])
  surfCoM = CoM(prod.bdys[surf[1]:surf[2]])
  
  v = molCoM - surfCoM |> (x -> repeat(x, N))  |> (x -> x / norm(x))

  p = dot(v, rc) |> abs |> (x -> round(x, digits=5))

  println(io, "SVP = $(p)")

  take!(io) |> String
end

