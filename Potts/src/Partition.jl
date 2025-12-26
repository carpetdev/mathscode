module Partition

using Polynomials
using LinearAlgebra: tr

function ising_part_periodic(n::Int)
    Ω = Iterators.product(Iterators.repeated((false, true), n)...)
    T₀ = zeros(Polynomial{BigInt}, 2^n, 2^n)
    # T₁ = zeros(typeof(x), 2^n, 2^n)

    for (i, Ωᵢ) ∈ zip(1:2^n, Ω),
        (j, Ωₒ) ∈ zip(1:2^n, Ω)

        # println(Ωᵢ)
        # println(Ωₒ)

        p = 0
        q = 0

        for r in 1:n
            p += Int(Ωᵢ[r] == Ωₒ[r])
            p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
            # q += Int(Ωₒ[r] == Ωₒ[mod1(r + 1, n)])
            # if r != n
            #     p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
            #     p += Int(Ωₒ[r] == Ωₒ[mod1(r + 1, n)])
            # end
        end

        # println(p)
        # q += p

        T₀[i, j] = Polynomial(1, p)
        # T₁[i, j] = x^q
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

    # return T₀
    partition = tr(T₀^n)
    # println(partition)
    # # partition = real_partition
    # # println(log(2, sum(partition)))
    # partition_roots = roots(partition)
    # display(scatter(partition_roots, aspect_ratio=:equal))

    return partition
end

function calculate_symmetry_classes(n::Int)
    Ω = Iterators.product(Iterators.repeated((false, true), n)...)
    classes::Vector{Vector{NTuple{n,Bool}}} = []
    seen_configs = Set{NTuple{n,Bool}}()
    for config in Ω
        if config in seen_configs
            continue
        end
        push!(classes, [config])
        push!(seen_configs, config)
        SD₂ₙ = Array{Function}(undef, 2, 2, n)
        SD₂ₙ[1, 1, :] = [t -> circshift(t, r - 1) for r in 1:n]
        SD₂ₙ[2, 1, :] = [t -> .!circshift(t, r - 1) for r in 1:n]
        SD₂ₙ[1, 2, :] = [t -> reverse(circshift(t, r - 1)) for r in 1:n]
        SD₂ₙ[2, 2, :] = [t -> .!reverse(circshift(t, r - 1)) for r in 1:n]

        for t in 1:2, s in 1:2, r in 1:n
            config′ = SD₂ₙ[t, s, r](config)
            if config′ in seen_configs
                continue
            end
            push!(classes[end], config′)
            push!(seen_configs, config′)
        end
    end
    return classes
end

function smul(S::Matrix{Polynomial{BigInt}}, T::Matrix{Polynomial{BigInt}}, n::Int)
    @assert length(symmetry_classes[n]), 2^n == size(S) == size(T)

    return
end

function spart(n)
    Ω = Iterators.product(Iterators.repeated((false, true), n)...)
    T₀ = zeros(Polynomial{Int,:x}, 2^n, 2^n)

    for (i, Ωᵢ) ∈ zip(1:2^n, Ω),
        (j, Ωₒ) ∈ zip(1:2^n, Ω)

        p = 0

        for r in 1:n
            p += Int(Ωᵢ[r] == Ωₒ[r])
            p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
        end

        T₀[i, j] = Polynomial(1, p)
    end

    partition = tr(T₀^n)

    return partition
end

end