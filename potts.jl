using Symbolics

@variables x

const n = 2
const T₀ = Matrix{typeof(x)}(undef, 2^n, 2^n)
const Ω = Iterators.product(Iterators.repeated((false, true), n)...)

for (i, Ωᵢ) ∈ zip(1:2^n, Ω),
    (j, Ωₒ) ∈ zip(1:2^n, Ω)

    p = 0

    for r ∈ 1:n
        p += Int(Ωᵢ[r] == Ωₒ[r])
        p += Int(Ωᵢ[r] == Ωᵢ[mod1(r, n)])
    end

    T₀[i, j] = x^p
end

display(T₀)