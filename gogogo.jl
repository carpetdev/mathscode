using Pkg

Pkg.activate("Potts")
Pkg.instantiate()

using Potts

for n in 5:12
    Save.part(n)
end