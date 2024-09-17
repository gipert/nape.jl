using Plots

# FIXME: this assumes main.jl has run

function plot_marginals()
    # TODO: please julia do not update x-axis range automatically
    plots = []
    for par in [:Γ12, :B]
        p = plot(samples, par, nbins=100)
        sup = support(prior[par])
        x = range(sup.lb, sup.ub, length=100)
        plot!(x, pdf(prior[par], x), color=:gray)
        push!(plots, p)
    end

    plot(plots..., size=(900, 400))
end

function plot_spectra()
    xlim = [1930, 2190]

    x = gerda.events.energy
    plot(
        x, fill(1E-3, length(x)), label="Data",
        st=:sticks, linewidth=2,
        markershape=:dtriangle, markerstrokewidth=0,
        xlim=xlim, ylim=(1E-4, 1E-1), yscale=:log10,
        xlabel="Energy [keV]", ylabel="BI [cts / (keV kg yr)]"
    )

    B_68 = BAT.smallest_credible_intervals(samples).B[1]

    plot!(
        xlim, fill(mean(B_68), 2),
        ribbon=width(B_68)/2,
        linewidth=0,
        fillalpha=0.2, label="BI"
    )

    plot!(
        xlim, fill(mode(samples).B, 2),
        st=:line, label=:none
    )

    P = partitions
    Γ12_90 = quantile(1E26 ./ samples.v.Γ12, 0.1)
    μs = sum(Γ12_90 * log(2) * N_A .* getfield.(P.ϵk, :val) / m_76 / ΔE)

    plot!(
        xlim, μs * pdf(Normal(2039, 2.5), xlim) .+ mode(samples).B,
        fill=:blue
    )
end
