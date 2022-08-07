# plot the errors
transD_GP.getchi2forall(opt)
ax = gcf().axes;
ax[2].set_ylim(15, 40)
χ² = aem.ndatalow + aem.ndatahigh
ax[2].plot(xlim(), [χ²/2 , χ²/2], "--", color="gray")
## plot the posterior resistivities
opt.xall[:] .= zall
transD_GP.plot_posterior(aem, opt, burninfrac=0.5, figsize=(5,6), nbins=50)
ax = gcf().axes
ax[1].invert_xaxis()
ax[1].step(log10.(ρ[2:end]), z[2:end], color="k", linewidth=3)
ax[1].step(log10.(ρ[2:end]), z[2:end], color="y", linewidth=1.5)
ax[1].set_ylim(280,0)
