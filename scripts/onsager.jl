using Polynomials
using Plots

N = 12
Plots.default(aspect_ratio=:equal, markersize=3, legend=false)

C(k, r) = cos(2π * k / N) + cos(2π * r / N)

onsageroots = Vector{ComplexF64}()
for k in 1:N, r in 1:N
    append!(onsageroots, roots(Polynomial([1, 2C(k, r), 2, -2C(k, r), 1])))
end

scatter(onsageroots)
