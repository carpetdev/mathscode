using Symbolics
using LinearAlgebra: tr
using PolynomialRoots
using Plots
import AMRVW

function potts(n)
    @variables x

    Ω = Iterators.product(Iterators.repeated((false, true), n)...)
    T₀ = zeros(typeof(x), 2^n, 2^n)
    T₁ = zeros(typeof(x), 2^n, 2^n)

    for (i, Ωᵢ) ∈ zip(1:2^n, Ω),
        (j, Ωₒ) ∈ zip(1:2^n, Ω)

        # println(Ωᵢ)
        # println(Ωₒ)

        p = 0
        q = 0 # amelia

        for r ∈ 1:n
            p += Int(Ωᵢ[r] == Ωₒ[r])
            p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
            q += Int(Ωₒ[r] == Ωₒ[mod1(r + 1, n)])
            # if r != n
            #     p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
            #     p += Int(Ωₒ[r] == Ωₒ[mod1(r + 1, n)])
            # end
        end

        # println(p)
        q += p

        T₀[i, j] = x^p
        T₁[i, j] = x^q
    end

    # configs = Iterators.product(Iterators.repeated((false, true), n^2)...)
    # real_partition = 0
    # for config ∈ configs
    #     p = 0
    #     config = reshape(collect(config), n, n)
    #     for i ∈ 1:n, j ∈ 1:n
    #         p += Int(config[i, j] == config[mod1(i + 1, n), j])
    #         p += Int(config[i, j] == config[i, mod1(j + 1, n)])
    #     end
    #     real_partition += x^p
    # end
    # println(real_partition)

    # println(sum(T₀))
    # display(T₀)

    partition = expand(tr(T₀^n))
    # println(partition)
    # partition = real_partition
    partition = Symbolics.coeff.(partition, x .^ collect(0:Symbolics.degree(partition)))

    partition_roots = AMRVW.roots(float.(partition))
    display(scatter(partition_roots))
end