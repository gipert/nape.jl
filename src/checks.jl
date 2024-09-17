function assert_in_fit_window(energy)
    E = energy
    @assert all(E .∈ [1930..2099] .|| E .∈ [2109..2114] .|| E .∈ [2124..2190])
end
