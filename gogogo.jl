using Pkg

Pkg.activate("Potts")
Pkg.instantiate()

using Potts

for n in 7:12
    Save.part(n)
end
