module Save

using FromFile
@from "Partition.jl" using Partition
@from "Load.jl" using Load

using JSON

function part(n)
    part = Partition.spart(n)
    JSON.json("data/parts/2_$(n)x$(n).json", part)
    return
end

function classes(n)
    symmetry = Symmetry{n}(Partition.calculate_symmetry_classes(n)...)
    JSON.json("data/symmetry/2_$(n)x$(n).json", symmetry)
    return
end

end