using Symbolics

@variables x

n = 2
T₀ = Matrix{typeof(x)}(undef, 2^n, 2^n)
Ω = Iterators.product(Iterators.repeated((false, true), n)...)

for (i, Ωᵢ) ∈ zip(1:2^n, Ω),
    (j, Ωₒ) ∈ zip(1:2^n, Ω),



end