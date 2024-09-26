using StatsBase: weights, Histogram
import PythonPlot as plt

# FIXME: this assumes main.jl has run! it's supposed to be run interactively
# after the fit

# this function calculates the observable counts density at a certain energy
# and given the samples. "density" means that there is no integral over energy
# (i.e. no multiplication by ΔE)
function counts_density_obs(partitions, energy, Γ12, B, α; bkg=true, signal=true)
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

    μk = (bkg ? μbk : 0) .+ (signal ? μsk .* pdf.(Normal.(Qββ .+ _Δk, _σk), energy) : 0)

    return sum(μk)
end

# this function calculates the posterior of the counts density/(kg yr) at a certain energy
function counts_density_kgyr_l200_post(energy, samples; kwargs...)
    @info "calculating posterior for energy $energy"
    partitions = data[:legend200].partitions
    f = s -> counts_density_obs(partitions, energy, s.Γ12, s.legend200_B, s.α; kwargs...)
    return broadcast(f, samples.v) ./ sum(partitions.exposure)
end

# this function calculates the counts density/(kg yr) at the marginal mode
function counts_density_kgyr_l200_mode(energy, mode; kwargs...)
    partitions = data[:legend200].partitions
    m = mode
    n_cts = counts_density_obs(partitions, energy, m.Γ12, m.legend200_B, m.α; kwargs...)
    return n_cts ./ sum(partitions.exposure)
end

function plot_l200_result()
    # we don't really need all samples...
    _samples = samples[1:20_000]
    _weights = weights(_samples.weight)

    # for the background, we just need one (random) energy point
    # we'll compute the 68% smallest CI
    posterior = counts_density_kgyr_l200_post(2000, _samples, signal=false)
    _int = BAT.smallest_credible_intervals(posterior, _weights, nsigma_equivalent=1)
    length(_int) != 1 && @warn "[$E keV] 68% interval is disjoint"
    b_68 = first(_int)

    # let's use a narrow region around Qbb for the 90% CI on the signal
    energies = 2034:0.3:2043

    s_90 = Vector(undef, length(energies))

    # let's speed this up with multi-threading
    Threads.@threads for i in 1:length(energies)
        E = energies[i]
        posterior = counts_density_kgyr_l200_post(E, _samples, bkg=false)

        # no need for smalles interval here
        # _int = BAT.smallest_credible_intervals(posterior, _weights, nsigma_equivalent=1.64)
        s_90[i] = 0..quantile(posterior, _weights, 0.9)
    end

    # plot!

    fig, ax = plt.subplots(figsize=(5, 2.5))

    band = plt.fill_between(
        [1930, 2190], fill(b_68.left, 2), fill(b_68.right, 2),
        linewidth=0, alpha=0.5, color="#228833",
        label="68% C.I.",
        zorder=0,
    )

    b_mode = counts_density_kgyr_l200_mode(2000, mode(_samples), signal=false)
    line, = plt.plot(
        [1930, 2190], [b_mode, b_mode],
        color="#228833",
        zorder=1,
    )

    ax.set_xlim(1930, 2190)
    ax.set_xlabel("Energy [keV]", loc="right")
    ax.set_ylabel("Counts / (keV kg yr)")
    # ylim=(1E-4, 1E-1), yscale=:log10,

    limit_area = plt.fill_between(
        energies,
        getproperty.(s_90, :left), getproperty.(s_90, :right),
        linewidth=0, color="#4477AA",
        label="90% C.I.",
    )

    # plot the data
    unbinned = false

    # re-add that count in the gamma line
    x = [data[:legend200].events.energy..., 2118.465]
    data_art = nothing

    if unbinned
        data_art, stemlines, baseline = ax.stem(
            x, fill(0, length(x)), "#CC3311",
            basefmt="none",
        )
        data_art.set_markersize(1.5)
        # NOTE: comment the following if the stem has a length!
        data_art.set_zorder(99)
        data_art.set_clip_on(false)

        stemlines.set_linewidth(1)
        ax.set_ylim(0, 0.013)
    else
        w = 1/sum(data[:legend200].partitions.exposure)
        _, _, data_art = plt.hist(
            x,
            weights=fill(w, length(x)),
            bins=1930:1:2190,
            color="#CC3311",
        )
        ax.set_ylim(1E-4, 2)
        ax.set_yscale("log")
    end

    # now this mess to plot the excluded regions
    plt.hist(
        [2104, 2119],
        bins=[2099, 2109, 2114, 2124],
        weights=[9, 9],
        color="black", alpha=0.4,
    )

    ax.legend(
        [data_art, (band, line), limit_area],
        [
            raw"LEGEND-200 data [48.3 kg yr]",
            raw"$6.2^{+1.8}_{-2.5} \times 10^{-4}$ cts / (keV kg yr) [68% C.I.]",
            raw"$T^{0\nu}_{1/2} > 1.8 \times 10^{26}$ yr [90% C.I.]",
        ],
        loc="upper left",
        fontsize="small",
    )

    plt.tight_layout()
    fig.savefig("l200-result.pdf", bbox_inches="tight")

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
