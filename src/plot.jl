using Plots

# FIXME: this assumes main.jl has run, needs to be rewritten

function n_counts_perkev_obs(partitions::Table, energy, Γ12, B, α)::Float64
    Qββ = 2039.04 # keV
    N_A = 6.02214
    m_76 = 75.9214 # g/mol
    P = partitions

    ϵk = k -> getfield.(P.ϵk, k)

    μsk = Γ12 * log(2) * N_A * P.exposure .* (ϵk(:val) .+ α * ϵk(:err)) / m_76
    μbk = B * P.exposure

    Δk = k -> getfield.(P.Δk, k)
    σk = k -> getfield.(P.σk, k)

    # NOTE: this assumes that the posterior is equal to the prior
    _Δk = rand.(Normal.(Δk(:val), Δk(:err)))
    _σk = rand.(Normal.(σk(:val), σk(:err)))

    μk = μbk .+ μsk .* pdf.(Normal.(Qββ .+ _Δk, _σk), energy)

    return sum(μk)
end

function n_counts_l200_post(energy, samples)
    partitions = data[:legend200].partitions
    f = s -> n_counts_perkev_obs(partitions, energy, s.Γ12, s.legend200_B, s.α)
    return broadcast(f, samples.v) ./ sum(partitions.exposure)
end

function n_counts_l200_mode(energy, mode)
    partitions = data[:legend200].partitions
    m = mode
    n_cts = n_counts_perkev_obs(partitions, energy, m.Γ12, m.legend200_B, m.α)
    return n_cts ./ sum(partitions.exposure)
end

function plot_l200_result()

    energies = [1930, collect(2030:0.5:2043)..., 2190]
    i68 = Vector(undef, length(energies))
    i95 = Vector(undef, length(energies))

    Threads.@threads for i in 1:length(energies)
        @info "calculating posterior for energy $(energies[i])"
        posterior = n_counts_l200_post(energies[i], samples)
        i68[i] = first(BAT.smallest_credible_intervals(posterior))
        i95[i] = first(BAT.smallest_credible_intervals(posterior, nsigma_equivalent=2))
    end

    plot(
        xlim=[2020, 2060], ylim=(1E-4, 1E-1), yscale=:log10,
        xlabel="Energy [keV]", ylabel="BI [cts / (keV kg yr)]"
    )

    plot!(
        energies,
        getproperty.(i95, :left), fillrange=getproperty.(i95, :right),
        color="orange", fillalpha=0.3, linewidth=0,
    )

    plot!(
        energies,
        getproperty.(i68, :left), fillrange=getproperty.(i68, :right),
        color="green", fillalpha=0.3, linewidth=0,
    )

    plot!(
        energies,
        broadcast(e -> n_counts_l200_gmode(e, mode(samples)), energies),
        color="black"
    )

    x = data[:legend200].events.energy

    plot!(
        x, fill(2E-3, length(x)), label="Data",
        st=:sticks, linewidth=2,
        markershape=:dtriangle, markerstrokewidth=0,
        color="black",
    )
end

# function plot_marginals()
#     # TODO: please julia do not update x-axis range automatically
#     plots = []
#     for par in [:Γ12, :B]
#         p = plot(samples, par, nbins=100)
#         sup = support(prior[par])
#         x = range(sup.lb, sup.ub, length=100)
#         plot!(x, pdf(prior[par], x), color=:gray)
#         push!(plots, p)
#     end

#     plot(plots..., size=(900, 400))
# end

# function plot_spectra()
#     xlim = [1930, 2190]

#     x = gerda.events.energy
#     plot(
#         x, fill(1E-3, length(x)), label="Data",
#         st=:sticks, linewidth=2,
#         markershape=:dtriangle, markerstrokewidth=0,
#         xlim=xlim, ylim=(1E-4, 1E-1), yscale=:log10,
#         xlabel="Energy [keV]", ylabel="BI [cts / (keV kg yr)]"
#     )

#     B_68 = BAT.smallest_credible_intervals(samples).B[1]

#     plot!(
#         xlim, fill(mean(B_68), 2),
#         ribbon=width(B_68)/2,
#         linewidth=0,
#         fillalpha=0.2, label="BI"
#     )

#     plot!(
#         xlim, fill(mode(samples).B, 2),
#         st=:line, label=:none
#     )

#     P = partitions
#     Γ12_90 = quantile(1E26 ./ samples.v.Γ12, 0.1)
#     μs = sum(Γ12_90 * log(2) * N_A .* getfield.(P.ϵk, :val) / m_76 / ΔE)

#     plot!(
#         xlim, μs * pdf(Normal(2039, 2.5), xlim) .+ mode(samples).B,
#         fill=:blue
#     )
# end
