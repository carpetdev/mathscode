using Polynomials
using Plots
using LinearAlgebra: tr

# Alternative root finders
import PolynomialRoots
import AMRVW

default(aspect_ratio=:equal)

# struct TransferMatrix <: AbstractMatrix{Tuple{Int}}
#     matrix::Matrix{Tuple{Int}}

#     Base.zero(::Type{Tuple{Int}}) = (0,)

#     function TransferMatrix(n::Int)
#         new(zeros(Tuple{Int}, n, n))
#     end

#     Base.size(T::TransferMatrix) = 2
# end

function potts(n)
    Ω = Iterators.product(Iterators.repeated((false, true), n)...)
    T₀ = zeros(Polynomial{Int,:x}, 2^n, 2^n)
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

    partition = tr(T₀^n)
    # println(partition)
    # # partition = real_partition
    # # println(log(2, sum(partition)))
    # partition_roots = roots(partition)
    # display(scatter(partition_roots, aspect_ratio=:equal))

    return partition
end

function bairstow(P::Polynomial)
    roots = Vector{ComplexF64}()
    while (n = degree(P)) > 1
        # println("Poly ", P)
        u::BigFloat = P[n-1] != 0 ? P[n-1] / P[n] : 1 / P[n]
        v::BigFloat = P[n-2] != 0 ? P[n-2] / P[n] : 1 / P[n]
        step = [Inf, Inf]
        while sum(abs2, step) > 0.00001
            # println("NR ", u, ' ', v)
            b = Vector(undef, n + 1)
            b[begin+n] = b[begin+n-1] = 0
            f = Vector(undef, n + 1)
            f[begin+n] = f[begin+n-1] = 0
            for i in n-2:-1:0
                b[begin+i] = P[i+2] - u * b[begin+i+1] - v * b[begin+i+2]
                f[begin+i] = b[begin+i+2] - u * f[begin+i+1] - v * f[begin+i+2]
            end
            c = P[1] - u * b[begin+0] - v * b[begin+1]
            d = P[0] - v * b[begin+0]
            g = b[begin+1] - u * f[begin+0] - v * f[begin+1]
            h = b[begin+0] - v * f[begin+0]
            step = 1 / (v * g^2 + h * (h - u * g)) * [-h g; -g*v g*u-h] * [c, d]
            u, v = [u, v] - step
        end
        b = Polynomial([v, u, 1])
        append!(roots, AMRVW.roots(coeffs(b)))
        a = P
        Q = 0
        while (m = degree(a)) > 1
            q = Polynomial(a[m], m - 2)
            # println(m, ' ', degree(a))
            # println("Euclid progress ", a, ' ', q)
            Q += q
            a -= q * b
            chop!(a)
        end
        # println("Euclid error ", a)
        P = Q
    end
    append!(roots, AMRVW.roots(coeffs(P)))
    return roots
end