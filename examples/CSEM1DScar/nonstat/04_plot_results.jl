GeophysOperator.getchi2forall(opt, fsize=8, alpha=1)
ax = gcf().axes;
ax[2].set_xlim(2000001, 4000001)
ax[2].set_ylim(2050, 2310)
ax[2].ticklabel_format(style="plain")
plt.tight_layout()
gcf().text(0.02, 0.9, "a.", fontsize=14, color="red")
savefig("csem_scar_conv_ns_1.png", dpi=300)
GeophysOperator.getchi2forall(optlog10λ, fsize=8, alpha=1)
ax = gcf().axes;
ax[2].set_xlim(2000001, 4000001)
ax[2].set_ylim(2050, 2310)
ax[2].ticklabel_format(style="plain")
plt.tight_layout()
gcf().text(0.02, 0.9, "b.", fontsize=14, color="red")
savefig("csem_scar_conv_ns_2.png", dpi=300)
opt.xall[:] .= zall
GeophysOperator.plot_posterior(csem, opt, optlog10λ, burninfrac=0.5, fsize=8,
    figsize=(7.8,4), cmappdf="inferno", nbins=50, vmaxpc=0.6,
    plotmode=true, plotquantile=false)
ax = gcf().axes
ax[1].step(log10.(ρ[2:end]), z[2:end], alpha=0.4, color="w")
ax[1].axis([-0.5, 1.5, 1500, 2500])
ax[1].set_xticks(-0.5:0.5:1.5)
ax[1].invert_yaxis()
gcf().text(0.02, 0.9, "b.", fontsize=14, color="red")
savefig("csem_scar_post_ns.png", dpi=300)
M = GeophysOperator.assembleTat1(opt, :fstar, temperaturenum=1)
Random.seed!(10)
GeophysOperator.plotmodelfield!(csem, M[randperm(length(M))[1:50]],
                                       dz=dz, extendfrac=extendfrac, onesigma=false)
ax = gcf().axes
ax[1].set_ylim(0,3450)
ax[1].invert_yaxis()
gcf().text(0.02, 0.9, "b.", fontsize=14, color="red")
savefig("csem_scar_fwds_post_ns.png", dpi=300)
