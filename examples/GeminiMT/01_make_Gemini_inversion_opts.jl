using Revise, transD_GP, Distributed, DelimitedFiles

## properties that apply to both GPs

println("are we in business?")

# load the resistivity mesh centroids
xall = readdlm("examples/GeminiMT/resistivityMeshCentroids_Gemini.txt")
xall = xall'
ymin = 0.9*minimum(xall[1,:]); ymax = 0.9*maximum(xall[1,:]);
zmin = minimum(xall[2,:]); zmax = 0.9*maximum(xall[2,:]);
xbounds = [ymin ymax; zmin zmax]   #[y; z] (m)
pnorm = 2.

## make the (hollow) m2d operator

m2d_op =     transD_GP.m2d_op()

## make options for the multichannel lengthscale GP

nminlog10λ, nmaxlog10λ = 2, 250
fbounds                = [1.0 5.0; 1.0 5.0]  # (m), log units
δlog10λ = 0.3
λlog10λ = 0.1*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
sdev_proplog10λ = 0.05*diff(log10bounds, dims=2)[:]
sdev_poslog10λ = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
demean = false # essential for length scale change approximations
Klog10λ = transD_GP.GP.Mat32()

## Initialize a lengthscale model using these options
Random.seed!(12)
optlog10λ = transD_GP.OptionsStat(nmin = nminlog10λ,
                        nmax = nmaxlog10λ,
                        xbounds = xbounds,
                        fbounds = fbounds,
                        xall = xall,
                        λ = λlog10λ,
                        δ = δlog10λ,
                        demean = demean,
                        sdev_prop = sdev_proplog10λ,
                        sdev_pos = sdev_poslog10λ,
                        pnorm = pnorm,
                        quasimultid = false,
                        K = Klog10λ,
                        timesλ = 3.6
                        )

## make options for the nonstationary actual properties GP

nmin, nmax = 2, 500
log10ρbounds = [-1. 1.]     # ohm-m (log units)
demean_ns = true
sdev_prop = 0.05*diff(fbounds, dims=2)[:]
sdev_pos = 0.05*abs.([diff([extrema(xall[1,:])...])[1], diff([extrema(xall[2,:])...])[1]])
δ = 0.25
K = transD_GP.GP.Mat32()

## Initialize model for the nonstationary properties GP
Random.seed!(13)
opt = transD_GP.OptionsNonstat(optlog10λ,
                        nmin = nmin,
                        nmax = nmax,
                        fbounds = log10bounds,
                        δ = δ,
                        demean = demean_ns,
                        sdev_prop = sdev_prop,
                        sdev_pos = sdev_pos,
                        pnorm = pnorm,
                        K = K
                        )
