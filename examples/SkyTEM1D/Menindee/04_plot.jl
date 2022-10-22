dr = 20
idx = 1:100:length(soundings)
burninfrac = 0.5
vmin, vmax = -2.5, 0.5
transD_GP.summaryAEMimages(soundings, opt,
                   zstart=zstart,
                   dz=dz,
                   dr=dr,
                   cmap="jet",
                   extendfrac=extendfrac,
                   vmin=vmin, vmax=vmax,
                   idx=idx,
                   nlayers=nlayers,
                   useML=useML,
                   preferEright=true,
                   saveplot=true
                   )
##
computeforwards = true
nforwards = 10
nbins = 50
transD_GP.plotindividualAEMsoundings(soundings, aem, opt, idx;
    burninfrac, nbins, computeforwards,
    nforwards)

