module Potts

export Partition, Roots, Save, Load

using FromFile
@from "Partition.jl" using Partition
@from "Roots.jl" using Roots
@from "Save.jl" using Save
@from "Load.jl" using Load

# struct TransferMatrix <: AbstractMatrix{Tuple{Int}}
#     matrix::Matrix{Tuple{Int}}

#     Base.zero(::Type{Tuple{Int}}) = (0,)

#     function TransferMatrix(n::Int)
#         new(zeros(Tuple{Int}, n, n))
#     end

#     Base.size(T::TransferMatrix) = 2
# end

end