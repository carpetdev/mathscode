n = 6
t = (false, true, false, false, true, false)

SD = Array{Function}(undef, 2, 2, n)
SD[1, 1, :] = [t -> circshift(t, r - 1) for r in 1:n]
SD[2, 1, :] = [t -> .!circshift(t, r - 1) for r in 1:n]
SD[1, 2, :] = [t -> circshift(reverse(t), r - 1) for r in 1:n]
SD[2, 2, :] = [t -> .!circshift(reverse(t), r - 1) for r in 1:n]

i = 1
j = 1
for k in 1:n
    @show t |> SD[i, j, k] |> SD[i, j, k]
end