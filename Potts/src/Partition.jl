module Partition

using FromFile
@from "Load.jl" using Load

using Polynomials
using LinearAlgebra: tr
using Bijections
using OrderedCollections
using BitIntegers

const Polynomial = Polynomial{UInt256} # `Polynomials.SparsePolynomial(c,d)` and `Polynomials.SparseVectorPolynomials(c,d)` don't work (the first stack overflows and the second assumes `d=1`). The first can be fixed with `Polynomials.SparsePolynomial([c],d)`.

function part(q::Int, n::Int) # ising_part_periodic
    if q == 2
        Ω = Load.symmetry_class(n).ordered_configs
    else
        Ω = Iterators.product(Iterators.repeated(1:q, n)...)
    end
    T₀ = zeros(Polynomial, q^n, q^n)
    # T₁ = zeros(typeof(x), 2^n, 2^n)

    for (i, Ωᵢ) in zip(1:q^n, Ω),
        (j, Ωₒ) in zip(1:q^n, Ω)

        # println(Ωᵢ)
        # println(Ωₒ)

        p = 0
        # s = 0

        for r in 1:n
            p += Int(Ωᵢ[r] == Ωₒ[r])
            p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
            # s += Int(Ωₒ[r] == Ωₒ[mod1(r + 1, n)])
            # if r != n
            #     p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
            #     p += Int(Ωₒ[r] == Ωₒ[mod1(r + 1, n)])
            # end
        end

        # println(p)
        # s += p

        T₀[i, j] = Polynomial([1], p)
        # T₁[i, j] = x^s
    end

    partition = tr(T₀^n)
    # println(partition)
    # partition_roots = roots(partition)
    # display(scatter(partition_roots, aspect_ratio=:equal))

    return partition
end

function dihedral(n::Int)
    SD = Array{Function}(undef, 2, 2, n)
    SD[1, 1, :] = [t -> circshift(t, r - 1) for r in 1:n]
    SD[2, 1, :] = [t -> .!circshift(t, r - 1) for r in 1:n]
    SD[1, 2, :] = [t -> reverse(circshift(t, r - 1)) for r in 1:n]
    SD[2, 2, :] = [t -> .!reverse(circshift(t, r - 1)) for r in 1:n]
    return SD
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

        SD = dihedral(n)
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
    T = zeros(Polynomial, length(classes), 2^n)

    SD = dihedral(n)
    function invert((f, s, r)::NTuple{3,Int})
        if s == 1
            return (f, s, mod1(2 - r, n))
        elseif s == 2
            return (f, s, r)
        end
    end

    for (i, Ωᵢ) in zip(1:length(classes), reps),
        (j, sigma) in zip(1:2^n, Iterators.flatten(classes))

        Ωₒ = SD[sigma...](reps[class_enum[j][1]])
        p = 0
        for r in 1:n
            p += Int(Ωᵢ[r] == Ωₒ[r])
            p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
        end
        T[i, j] = Polynomial([1], p)
    end
    # return T

    m = n
    if m & 1 == 1
        Tⁿ = T
    else
        Tⁿ = zeros(Polynomial, length(classes), 2^n)
        rep_index = 1
        for (i, class) in enumerate(classes)
            Tⁿ[i, rep_index] = 1
            rep_index += length(class)
        end
    end

    function smul(S::Matrix{Polynomial}, T::Matrix{Polynomial}, n::Int)
        @assert (length(classes), 2^n) == size(S) == size(T)
        out = zeros(Polynomial, length(classes), 2^n) # What if we don't specify type
        Threads.@threads for (i, j) in collect(Iterators.product(1:length(classes), 1:2^n))
            for k in 1:2^n
                c, d = class_enum[k]
                a = config_by_index(SD[invert(classes[c][d])...](config_by_index[j]))
                out[i, j] += S[i, k] * T[c, a]
            end
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

function spart′(n::Int)
    (; classes, reps) = Load.symmetry_class(n)
    class_enum = [(c, d) for c in 1:length(classes) for d in 1:length(classes[c])]
    T = zeros(Polynomial, length(classes), 2^n)

    SD = dihedral(n)
    function invert((f, s, r)::NTuple{3,Int})
        if s == 1
            return (f, s, mod1(2 - r, n))
        elseif s == 2
            return (f, s, r)
        end
    end

    for (i, Ωᵢ) in zip(1:length(classes), reps),
        (j, sigma) in zip(1:2^n, Iterators.flatten(classes))

        Ωₒ = SD[sigma...](reps[class_enum[j][1]])
        p = 0
        for r in 1:n
            p += Int(Ωᵢ[r] == Ωₒ[r])
            p += Int(Ωᵢ[r] == Ωᵢ[mod1(r + 1, n)])
        end
        T[i, j] = Polynomial([1], p)
    end

    R = sum(T, dims=2)
    # display(Polynomials.Polynomial{BigInt}.(T))

    #region matrix-vector multiplication manual
    for i in 1:n-2
        # println(i)
        R_new = zeros(Polynomial, length(classes))
        Threads.@threads for i in 1:length(classes)
            for k in 1:2^n
                c, _ = class_enum[k]
                R_new[i] += T[i, k] * R[c] # Try just making R a vector with repeats so you can do matrix-vector mulitiplication
            end
        end
        R = R_new
    end
    partition = sum(R .* (length(class) for class in classes))
    #endregion

    #region matrix-vector multiplication alternative
    # for i in 1:n-2
    #     print(i)
    #     R = reduce(vcat, fill.(R, (length(class) for class in classes)))
    #     R = T * R
    # end
    # partition = sum(R)
    #endregion

    return partition
end

end