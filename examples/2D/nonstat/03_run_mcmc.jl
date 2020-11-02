## set up McMC
nsamples, nchains, nchainsatone = 12001, 4, 1
Tmax = 2.50
addprocs(nchains)
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using transD_GP
## run McMC
@time transD_GP.main(optlog10λ, opt, img, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)
rmprocs(workers())
