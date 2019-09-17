import MPI

# always initiate MPI with this:
MPI.Init()

comm = MPI.COMM_WORLD           # the default MPI communicator, which has all processors invoked with mpirun
rank = MPI.Comm_rank(comm)     # the id (index) of this processor in communicator comm
nprocs = MPI.Comm_size(comm)     # the total number of processors in the communicator group comm

# define number of PT chains, the temperature ladder, and the number of total MCMC samples
nsamples, nchains, Tmax = 4001, 2, 2.5
@assert mod(nprocs,nchains) == 0 # there should be the same whole number of procs per PT chain
# create a team for each PT chain
nProcPerTeam = Int(nprocs/nchains)

# this process' team number is:
team = Integer(floor(rank/nProcPerTeam))

# define the worker captains (the "leader" threads of each team)
WorkerCaptains = 0:nProcPerTeam:nprocs-1

if in(rank,WorkerCaptains)
    println("my rank is $rank, so I'm a worker captain!")
    any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
    using ImageRegression, TransD_GP, MCMC_Driver, Distributed
end

# let's wrap things up
MPI.Finalize()
