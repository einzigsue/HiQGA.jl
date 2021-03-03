
using MPIClusterManagers, Distributed

## set up McMC

# Initialize MPI and get communicator, rank and number of MPI processes in communicator:
#MPI.Init()
manager = MPIManager(np=4)
addprocs(manager)
println("Added procs $(procs())")
# make sure nproc (below) matches np=? (above)
@everywhere nsamples, nchains, nchainsatone, nproc = 6, 2, 1, 4
@assert mod(nproc,nchains) == 0 # there should be a whole number of procs per PT chain

Tmax = 1.2

# make sure everyone has the right modules
@everywhere using Distributed
@everywhere using transD_GP
@everywhere import MPI
@everywhere import Mare2dem
@everywhere using SparseArrays
@everywhere using SuiteSparse
@everywhere using Random
@everywhere using Printf

@everywhere include("src/m2d.jl")
@everywhere filenameroot = "prism1.0"
@everywhere nprocperchain = Int(nproc/nchains)
@mpi_do manager load_m2d(filenameroot,nprocperchain)
println("we made it out of the load_m2d call!")

@time transD_GP.main(optlog10λ, opt, m2d_op,
                     Tmax=Tmax, nsamples=nsamples, nchains=nchains, nchainsatone=nchainsatone,
                     m2d_flag=true)

rmprocs(workers())
println("Exiting")
exit()
