module Partition

using FromFile
@from "Load.jl" using Load

using Polynomials
using LinearAlgebra: tr
using Bijections

function ising_part_periodic(n::Int)
    Ω = Iterators.product(Iterators.repeated((false, true), n)...)
    T₀ = zeros(Polynomial{BigInt}, 2^n, 2^n)
    # T₁ = zeros(typeof(x), 2^n, 2^n)

    for (i, Ωᵢ) in zip(1:2^n, Ω),
        (j, Ωₒ) in zip(1:2^n, Ω)

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
    # for config in configs
    #     p = 0
    #     config = reshape(collect(config), n, n)
    #     for i in 1:n, j in 1:n
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
    classes = Vector{Vector{NTuple{2,Function}}}() # List of lists of a tuple of symmetry and inverse
    # reps = Vector{NTuple{n,Bool}}()
    seen_configs = Set{NTuple{n,Bool}}()
    config_by_index = Bijection{Int,NTuple{n,Bool}}()
    index = 1

    for config in Ω
        if config in seen_configs
            continue
        end
        push!(classes, [(identity, identity)])
        # push!(reps, config)
        push!(seen_configs, config)
        config_by_index[index] = config
        index += 1
        SD = Array{Function}(undef, 2, 2, n)
        SD[1, 1, :] = [t -> circshift(t, r - 1) for r in 1:n]
        SD[2, 1, :] = [t -> .!circshift(t, r - 1) for r in 1:n]
        SD[1, 2, :] = [t -> reverse(circshift(t, r - 1)) for r in 1:n]
        SD[2, 2, :] = [t -> .!reverse(circshift(t, r - 1)) for r in 1:n]

        for f in 1:2, s in 1:2, r in 1:n
            config′ = SD[f, s, r](config)
            if config′ in seen_configs
                continue
            end
            push!(classes[end], (SD[f, s, r], SD[f, s, mod1(2 - r, n)]))
            push!(seen_configs, config′)
            config_by_index[index] = config′
            index += 1
        end
    end
    return classes, config_by_index
    # return classes, reps, config_by_index
end

function smul(S::Matrix{Polynomial{BigInt}}, T::Matrix{Polynomial{BigInt}}, n::Int)
    # classes, reps, config_by_index = Load.symmetry_classes[n]
    classes, config_by_index = calculate_symmetry_classes(n)
    class_enum = [(i, j) for i in 1:length(classes) for j in 1:length(classes[i])]
    @assert (length(classes), 2^n) == size(S) == size(T)
    out = zeros(Polynomial{BigInt}, length(classes), 2^n)
    for i in 1:length(classes), j in 1:2^n, k in 1:2^n
        out[i, j] += S[i, k] * T[rep_index(k), config_by_index(classes[2]())]
    end
    return
end

function spart(n)
    Ω = Iterators.product(Iterators.repeated((false, true), n)...)
    T₀ = zeros(Polynomial{BigInt}, 2^n, 2^n)

    for (i, Ωᵢ) in zip(1:2^n, Ω),
        (j, Ωₒ) in zip(1:2^n, Ω)

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