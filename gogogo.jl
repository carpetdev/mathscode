using Pkg

Pkg.activate("Potts")
Pkg.instantiate()

using Potts

for n in 12:16
    Save.part(n)
end
