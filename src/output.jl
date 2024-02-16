
function printSVP(reac, ts, prod)
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

  take!(io) |> String
end

