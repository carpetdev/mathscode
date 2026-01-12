module Save

using FromFile
@from "Partition.jl" using Partition
@from "Load.jl" using Load

using JSON

function parts(n)
    parts = Load.parts
    part = Partition.ising_part_periodic(n)
    if n > length(parts) + 1
        JSON.json("data/parts/2_$(n)x$(n).json", part)
        return
    elseif n == length(parts) + 1
        push!(parts, part)
    else
        parts[n] = part
    end
    JSON.json("data/parts/2_nxn.json", parts)
    return
end

function classes(n)
    classes, reps, config_by_index = Partition.calculate_symmetry_classes(n)
    return
end

end