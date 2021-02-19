## set up McMC
using Distributed
nsamples, nchains, nchainsatone = 100001, 8, 1
Tmax = 3
addprocs(nchains)
##init packages on workers
@info "workers are $(workers())"
@everywhere using Distributed
@everywhere using transD_GP
## run McMC
@time transD_GP.main(opt, optdummy, optn, tempest, Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone)

## clean up
rmprocs(workers())

# ## test a move for misfit change
# z = [-100000.0, 0.0, 2.0, 4.120000000000002, 6.367199999999997, 8.749232000000003, 11.274185919999999, 13.9506370752, 16.787675299712006, 19.794935817694725, 22.982631966756408, 26.361589884761795, 29.943285277847504, 33.73988239451835, 37.76427533818946, 42.030131858480836, 46.55193976998968, 51.34505615618907, 56.4257595255604, 61.81130509709404, 67.51998340291968, 73.57118240709487, 79.98545335152056, 86.7845805526118, 93.99165538576851, 101.63115470891464, 109.72902399144951, 118.31276543093648, 127.41153135679268, 137.05622323820023, 147.27959663249226, 158.1163724304418, 169.60335477626833, 181.77955606284445, 194.68632942661512, 208.36750919221203, 222.8695597437448, 238.24173332836943, 254.53623732807162, 271.80841156775597, 290.1169162618213];
#
# ρ = [1.0e12, 11.518249922739678, 11.510479873302128, 11.4667122672664, 11.208562758178864, 9.75070498259659, 4.441154779564496, 0.5310449162001603, 69.83096849995604, 51.18267665590463, 13.755473422060344, 11.675713763440733, 11.528808707547665, 11.520495076484165, 11.520129613622801, 11.520117042532757, 11.5201167039629, 11.520116696859642, 11.520116696744546, 11.520116696743122, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311, 11.52011669674311]
#
# geovec = [-120.0, -76.0387905326985, -114.85633859374474, 0.0, 0.0, 1.8822087446787323, 0.0, 0.0, 0.0, 0.0]
#
# gv2 = [-120.0, -76.5, -114.85633859374474, 0.0, 0.0, 1.8822087446787323, 0.0, 0.0, 0.0, 0.0]
