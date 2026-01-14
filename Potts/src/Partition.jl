module Partition

using FromFile
@from "Load.jl" using Load

using Polynomials
using LinearAlgebra: tr
using Bijections
using OrderedCollections

function part(n::Int) # ising_part_periodic
    # Ω = Iterators.product(Iterators.repeated((false, true), n)...)
    Ω = Load.symmetry_class(n).ordered_configs
    T₀ = zeros(Polynomial{BigInt}, 2^n, 2^n)
    # T₁ = zeros(typeof(x), 2^n, 2^n)

    for (i, Ωᵢ) in zip(1:2^n, Ω),
        (j, Ωₒ) in zip(1:2^n, Ω)

        # println(Ωᵢ)
        # println(Ωₒ)

        p = 0
        # q = 0

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
    classes = Vector{Vector{NTuple{3,Int}}}() # List of lists of SD index for symmetry
    reps = Vector{NTuple{n,Bool}}()
    seen_configs = OrderedSet{NTuple{n,Bool}}()

    for config in Ω
        if config in seen_configs
            continue
        end
        push!(classes, [((1, 1, 1))])
        push!(reps, config)
        push!(seen_configs, config)
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
            push!(classes[end], (f, s, r))
            push!(seen_configs, config′)
        end
    end
    @assert sum(length(class) for class in classes) == length(seen_configs) == 2^n
    return classes, reps, collect(seen_configs)
end

function spart(n::Int)
    (; classes, reps, ordered_configs) = Load.symmetry_class(n)
    class_enum = [(c, d) for c in 1:length(classes) for d in 1:length(classes[c])]
    config_by_index = Bijection([i => v for (i, v) in enumerate(ordered_configs)])
    T = zeros(Polynomial{BigInt}, length(classes), 2^n)
    SD = Array{Function}(undef, 2, 2, n)
    SD[1, 1, :] = [t -> circshift(t, r - 1) for r in 1:n]
    SD[2, 1, :] = [t -> .!circshift(t, r - 1) for r in 1:n]
    SD[1, 2, :] = [t -> reverse(circshift(t, r - 1)) for r in 1:n]
    SD[2, 2, :] = [t -> .!reverse(circshift(t, r - 1)) for r in 1:n]

    for (i, Ωᵢ) in zip(1:length(classes), reps),
        (j, sigma) in zip(1:2^n, Iterators.flatten(classes))

        Ωₒ = SD[sigma...](reps[class_enum[j][1]])
        p = 0
        for r in 1:n
            p += Int(Ωᵢ[r] == Ωₒ[r])
            p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
        end
        T[i, j] = Polynomial(1, p)
    end
    # return T

    m = n
    if m & 1 == 1
        Tⁿ = T
    else
        Tⁿ = zeros(Polynomial{BigInt}, length(classes), 2^n)
        rep_index = 1
        for (i, class) in enumerate(classes)
            Tⁿ[i, rep_index] = 1
            rep_index += length(class)
        end
    end

    function smul(S::Matrix{Polynomial{BigInt}}, T::Matrix{Polynomial{BigInt}}, n::Int)
        class_enum = [(c, d) for c in 1:length(classes) for d in 1:length(classes[c])]
        @assert (length(classes), 2^n) == size(S) == size(T)
        out = zeros(Polynomial{BigInt}, length(classes), 2^n)
        for i in 1:length(classes), j in 1:2^n, k in 1:2^n
            c, d = class_enum[k]
            f, s, r = classes[c][d]
            out[i, j] += S[i, k] * T[c, config_by_index(SD[f, s, mod1(2 - r, n)](config_by_index[j]))]
        end
        return out
    end
    while m > 0
        m >>= 1
        T = smul(T, T, n)
        if m & 1 == 1
            Tⁿ = smul(Tⁿ, T, n)
        end
    end

    partition = 0
    rep_index = 1
    for (i, class) in enumerate(classes)
        partition += length(class) * Tⁿ[i, rep_index]
        rep_index += length(class)
    end
    return partition
end

end