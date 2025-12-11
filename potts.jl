using Polynomials
using Plots
using LinearAlgebra: tr
using Serialization
using ProgressLogging
using Random
using Attractors
using CairoMakie

# Alternative root finders
import PolynomialRoots
import AMRVW

default(aspect_ratio=:equal, markersize=3, legend=false)

const parts = deserialize("parts.dat")

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
    P₀ = P
    preroots = Vector{ComplexF64}()
    while (n = degree(P)) > 1
        # println("Poly ", P)
        u::BigFloat = P[n-1] != 0 ? P[n-1] / P[n] : 1 / P[n]
        v::BigFloat = P[n-2] != 0 ? P[n-2] / P[n] : 1 / P[n]
        step = [Inf, Inf]
        while sum(abs2, step) > 1e-25
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
        root = quadraticroot(Polynomial([v, u, 1]))
        # root = newtonraphson(P₀, derivative(P₀), root)
        push!(preroots, root)
        a = P₀
        b = fromroots(preroots)
        b *= conj(b)
        @show b
        Q = 0
        while degree(a) >= degree(b)
            @show degree(a)
            q = Polynomial(a[degree(a)], degree(a) - degree(b))
            # println(m, ' ', degree(a))
            # println("Euclid progress ", a, ' ', q)
            Q += q
            a -= q * b
            chop!(a)
        end
        println("Euclid error ", a)
        @show Q
        P = Q
        for i in 0:degree(P)
            if imag(P[i]) > 0.1
                println("HERE", imag(P[i]))
            end
            P[i] = real(P[i])
        end
    end
    println("Final Poly: ", P)
    r = Polynomials.roots(P)
    if !isempty(r)
        push!(preroots, r[1])
    end

    roots = Vector{ComplexF64}()
    colours = Vector{Symbol}()
    for root in preroots
        # root = newtonraphson(P₀, derivative(P₀), root)
        error = P₀(root)
        if abs(root) < 0.1
            println(root)
        end
        if abs(error) > 1e-10
            # println(error)
            push!(colours, :red, :red)
        else
            push!(colours, :black, :black)
        end
        push!(roots, root, conj(root))
    end
    @show degree(P₀) length(roots) length(unique(roots))

    display(scatter(roots, seriescolor=colours))
    return roots
end

function newtonraphson(f, f′, initial; tol=1e-16, maxIter=1e5)
    x = Complex{BigFloat}(initial)
    fx = f(x)
    iter = 0
    while abs(fx) > tol && iter < maxIter
        x -= fx / f′(x)
        fx = f(x)
        iter += 1
    end
    return x
end

function quadraticroot(P::Polynomial)
    a, b, c = P[2], P[1], P[0]
    return -b / 2a + √Complex((b / 2a)^2 - c / a)
end

function sub(P::Polynomial, s::Number)
    P = convert(Polynomial{typeof(s)}, copy(P))
    for i in 0:degree(P)
        P[i] *= s^i
    end
    return P
end

function hubbard(P::Polynomial; ε=1e-10, R=3, rootcount=degree(P))
    roots = Vector{Complex{BigFloat}}()
    debug_roots = []
    nonconv = Vector{Complex{BigFloat}}()
    debug_nonconv = []
    debug_dupe = []
    haszeroroot = false
    if P[0] == 0
        P = Polynomial(P[findfirst(!isequal(0), P):end])
        push!(roots, 0)
        haszeroroot = true
    end
    d = degree(P)
    s = ceil(Int, 0.26632log(d))
    N = ceil(Int, 8.32547d * log(d))
    @show d s N
    # R = 1 + maximum(abs(P[i]) / abs(P[d]) for i in 0:d-1)

    # max_N′ = N ÷ 2 + 1 # Just upper half plane - consider other symmetries
    max_N′ = N
    N′ = zeros(Int, max_N′) # More interesting enumeration of N
    i = 1
    for r in 0:d-1, q in 0:(max_N′-1-r)÷d
        N′[i] = q * d + r
        i += 1
    end

    S = (R * (1 + √(big(2))) * (big(d - 1) / d)^(big(2v - 1) / 4s) * exp(im * big(π) * big(2j - v % 2) / N) for j in N′, v in 1:s, R in [0.8, 3])
    P′ = derivative(P)
    K = ceil(Int, d * log(R / ε))
    # K = 10_000
    @show K
    @progress for point in S
        if length(roots) >= rootcount
            return roots
        end
        if length(roots) >= d + Int(haszeroroot) # Is Z square-free?
            break
        end
        point₀ = point
        hasconverged = false
        for _ in 1:K
            update = P(point) / P′(point)
            point -= update
            if abs(update) < ε / d
                hasconverged = true
                break
            end
        end
        if hasconverged
            isnewroot = true
            for root in roots
                if abs(root - point) < ε
                    isnewroot = false
                    break
                end
            end
            if isnewroot
                push!(roots, point, conj(point))
                push!(debug_roots, Dict(:initial => point₀, :final => point))

                @show point ComplexF64(P(point)) length(roots)
            else
                push!(debug_dupe, Dict(:initial => point₀, :final => point))
            end
        else
            push!(nonconv, point)
            push!(debug_nonconv, Dict(:initial => point₀, :final => point))
        end
    end
    display(scatter(roots))
    @show d length(roots) + length(nonconv)
    @show length(debug_dupe)
    serialize("debug_roots.dat", debug_roots)
    serialize("debug_dupe.dat", debug_dupe)
    serialize("debug_nonconv.dat", debug_nonconv)
    return roots, nonconv, debug_roots, debug_dupe, debug_nonconv
end

function newton_map(z, p, n)
    z1 = z[1] + im * z[2]
    dz1 = newton_f(z1, p[1]) / newton_df(z1, p[1])
    z1 = z1 - dz1
    return SVector(real(z1), imag(z1))
end
part = copy(parts[6])
newton_f(x, p) = part(x)
newton_df(x, p) = derivative(part)(x)

ds = DiscreteDynamicalSystem(newton_map, [0.1, 0.2], [3.0])
xg = yg = range(-2.5, 2.5; length=400)
grid = (xg, yg)
# Use non-sparse for using `basins_of_attraction`
mapper_newton = AttractorsViaRecurrences(ds, grid;
    sparse=false, consecutive_lost_steps=1000
)
basins, attractors = basins_of_attraction(mapper_newton; show_progress=false)
heatmap_basins_attractors(grid, basins, attractors; markers=Dict(i => :circle for i in 1:length(attractors)))