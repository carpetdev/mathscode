module Load

export Symmetry

using Polynomials
using JSON
using Bijections

struct Symmetry{n}
    classes::Vector{Vector{NTuple{3,Int}}}
    reps::Vector{NTuple{n,Bool}} # Can also just find this from ordered_configs and classes together
    ordered_configs::Vector{NTuple{n,Bool}}
end

function part(n::Int)
    return JSON.parsefile("data/parts/2_$(n)x$(n).json", Polynomial{BigInt})
end

function part′(n::Int)
    return JSON.parsefile("data/parts/2_$(n)x$(n)'.json", Polynomial{BigInt})
end

function symmetry_class(n::Int)
    return JSON.parsefile("data/symmetry/2_$(n)x$(n).json", Symmetry{n})
end

end