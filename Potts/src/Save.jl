module Save

using FromFile
@from "Partition.jl" include Partition

using JSON3
using JLD2

function parts(n)
    parts = Partition.parts
    part = Partition.ising_part_periodic(n)
    @show part
    if n > length(parts)
        JSON3.write("data/parts/2_$(n)x$(n).json", part)
        return
    elseif n == length(parts)
        push!(parts, part)
    else
        parts[n] = part
    end
    JSON3.write("data/parts/2_nxn.json", parts)
    return
end

function classes(n)
    classes = Partition.calculate_symmetry_classes(n)
    @save "data/symmetry/$(n)x$(n).jld2" classes
end

end