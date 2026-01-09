module Save

using FromFile
@from "Partition.jl" using Partition
@from "Load.jl" using Load

using JSON3
using JLD2

function write(path::String, list::Vector{T}, item::T, n::Int, format::Symbol=:json) where T
    if n > length(list) + 1
        JSON3.write(path * "2_$(n)x$(n).json", item)
        return
    elseif n == length(list) + 1
        push!(list, item)
    else
        list[n] = item
    end
    JSON3.write(path * "2_nxn.json", list)


end

function parts(n)
    parts = Load.parts
    part = Partition.ising_part_periodic(n)
    if n > length(parts) + 1
        JSON3.write("data/parts/2_$(n)x$(n).json", part)
        return
    elseif n == length(parts) + 1
        push!(parts, part)
    else
        parts[n] = part
    end
    JSON3.write("data/parts/2_nxn.json", parts)
    return
end

function classes(n) # remove reps?
    classes, reps, config_by_index = Partition.calculate_symmetry_classes(n)
    @save "data/symmetry/$(n)x$(n).jld2" (classes, reps, config_by_index)
end

end