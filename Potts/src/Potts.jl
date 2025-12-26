module Potts

export Partition, Roots, Save, Load

using FromFile
@from "Partition.jl" import Partition
@from "Roots.jl" import Roots
@from "Save.jl" import Save
@from "Load.jl" import Load

# struct TransferMatrix <: AbstractMatrix{Tuple{Int}}
#     matrix::Matrix{Tuple{Int}}

#     Base.zero(::Type{Tuple{Int}}) = (0,)

#     function TransferMatrix(n::Int)
#         new(zeros(Tuple{Int}, n, n))
#     end

#     Base.size(T::TransferMatrix) = 2
# end

end