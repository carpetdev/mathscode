module Load

using Polynomials
using JSON
using Bijections

struct Part
    n::Int
    polynomial::Polynomial{BigInt}
end

struct Symmetry
    classes::Vector{Vector{NTuple{2,NTuple{3,Int}}}}
    reps::Vector{NTuple{n,Bool}}
    config_by_index::Bijection{Int,NTuple{n,Bool}}
end

const parts = JSON.parsefile("data/parts/2_nxn.json", Vector{Part})

end