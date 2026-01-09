module Load

using Polynomials
using JSON3

const parts = JSON3.read("data/parts/2_nxn.json", Vector{Polynomial{BigInt}})

end