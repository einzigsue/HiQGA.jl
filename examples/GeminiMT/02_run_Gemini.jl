
using MPIClusterManagers, Distributed

## set up McMC
nsamples, nchains, nchainsatone = 6, 2, 1
Tmax = 1.2

# Initialize MPI and get communicator, rank and number of MPI processes in communicator:
#MPI.Init()
manager = MPIManager(np=(nchains*2))
addprocs(manager)
println("Added procs $(procs())")

# make sure everyone has the right modules
@everywhere using Distributed
@everywhere using transD_GP
@everywhere import MPI
@everywhere import Mare2dem
@everywhere using SparseArrays
@everywhere using SuiteSparse
@everywhere using Random
@everywhere using Printf

@time transD_GP.main(optlog10λ, opt, m2d_op,
                     Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone,
                     m2d_flag=true)

rmprocs(workers())
println("Exiting")
exit()
