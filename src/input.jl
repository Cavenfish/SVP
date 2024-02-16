abstract type Particle end 

struct Atom <: Particle
  s::Char
  m::Float64
  r::Vector{Float64}
end

struct Hess
  hessian::Matrix{Float64}
  modes::Matrix{Float64}
  freq::Vector{Float64}
  bdys::Vector{Particle}
end

function readAtoms(data; ext=".hess")
  n = popfirst!(data) |> (x -> parse(Int, x))
  
  if ext == ".xyz"
    popfirst!(data)
  end

  bdys = Vector{Particle}(undef, n)
  
  for i = 1:n
    tmp = popfirst!(data) |> split
    s   = tmp[1][1]
    m   = parse(Float64, tmp[2])
    r   = parse.(Float64, tmp[3:5])
    
    bdys[i] = Atom(s, m, r)
  end

  bdys
end

function readVector(data)
  s   = popfirst!(data) |> split |> (x -> parse(Int, x[1]))
  vec = zeros(s)

  for i in 1:s
    vec[i] = popfirst!(data) |> split |> (x -> parse(Float64, x[2]))
  end

  vec
end

function readMatrix(data)
  s   = popfirst!(data) |> split |> (x -> parse.(Int, x))
  length(s) == 2 ? (l,m) = s : (l,m) = repeat(s, 2) 
  mat = zeros(l,m)

  J = popfirst!(data) |> split |> (x -> parse.(Int, x) .+ 1)
  
  while !(isempty(J))
    for i in 1:m
      tmp = popfirst!(data) |> split |> (x -> parse.(Float64, x[2:end]))

      mat[i, J] .= tmp
    end
    J = popfirst!(data) |> split |> (x -> parse.(Int, x) .+ 1)
  end
  
  mat
end

function readHess(file::String)
  data  = readlines(file)

  i     = findfirst(e -> occursin("\$hessian", e), data)
  hessi = readMatrix(data[i+1:end])

  i     = findfirst(e -> occursin("normal_modes", e), data)
  modes = readMatrix(data[i+1:end])

  i     = findfirst(e -> occursin("\$vibrational_frequencies", e), data)
  freq  = readVector(data[i+1:end])

  i     = findfirst(e -> occursin("\$atoms", e), data)
  bdys  = readAtoms(data[i+1:end])

  Hess(hessi, modes, freq, bdys)
end