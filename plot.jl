using Plots

function plot_marginals()
    plot(plot(samples, :Γ12, nbins=100), plot(samples, :B, nbins=100), size=(900, 400))
end

function plot_spectra()
    xlim = [1930, 2160]

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
        fillalpha=0.2, label=:none
    )

    # plot!(
    #     xlim, fill(mode(samples).B, 2),
    #     st=:line
    # )
end
