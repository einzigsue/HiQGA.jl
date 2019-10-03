import MPI
# make sure this directory is added to the LOAD_PATH
any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
# make sure every process has this module loaded
using MCMC_Driver

# always initiate MPI with this:
MPI.Init()

comm = MPI.COMM_WORLD           # the default MPI communicator, which has all processors invoked with mpirun
rank = MPI.Comm_rank(comm)     # the id (index) of this processor in communicator comm
nprocs = MPI.Comm_size(comm)     # the total number of processors in the communicator group comm

# define number of PT chains, the temperature ladder, and the number of total MCMC samples
nsamples, nchains, Tmax = 1, 4, 2.5
@assert mod(nprocs,nchains) == 0 # there should be the same whole number of procs per PT chain
# create a team for each PT chain
nProcPerTeam = Int(nprocs/nchains)

# this process' team number is:
team = Integer(floor(rank/nProcPerTeam))
# this process' team rank
teamRank = mod(rank,nProcPerTeam)

# define the worker captains (the "leader" threads of each team)
WorkerCaptains = Array{Int64,1}(0:nProcPerTeam:nprocs-1)
# if this process is a WorkerCaptain, find out which captain it is
if in(rank,WorkerCaptains)
    tmp = findall(x -> x == rank, WorkerCaptains)
    WCindex = tmp[1]
else
    WCindex = -1
end

# save the relevant MPI parameters (rank, etc) so every process remembers what's what
myMPIparams = MCMC_Driver.MPIparams(rank,nprocs,team,teamRank,WCindex,WorkerCaptains)
#=if(rank == 0)
    println("Total processors: $nprocs")
    println("Number of teams: $nchains")
    println("Processors per team: $nProcPerTeam")
println("my team: $team")
println("my team rank: $teamRank")=#

# each captain will load the data and metadata
if in(rank,WorkerCaptains)
    println("my rank is $rank, so I'm a worker captain!")

    # load the necessary modules
    using ImageRegression, TransD_GP

    # metadata about the problem
    img = ImageRegression.Img(
              filename         = "4.2.01.tiff",
              dx               = 10.0,
              fractrain        = 0.02,
              dec              = 2,
              gausskernelwidth = 7)
    # load the data
    d, sd, ftrain, Xtrain =  ImageRegression.get_training_data(img,
                       sdmaxfrac = 0.05,
                       ybreak = 1000,
                       takeevery = 4)
    # have one of the captains plot the data
    if rank == 0
        ImageRegression.plot_data(ftrain, Xtrain, img)
    end

    # define the points where we want to predict the image
    Xall = ImageRegression.get_all_prediction_points(img)

    # the MCMC algorithm options/parameters
    opt_in = TransD_GP.Options(
                  nmin = 2,
                  nmax = 20,
                  xbounds           = [img.x[1] img.x[end];img.y[1] img.y[end]],
                  fbounds           = [-1 2],
                  xall              = Xall,
                  λ                 = [150.0, 150.0],
                  δ                 = 0.2,
                  demean            = true,
                  save_freq         = 500,
                  dispstatstoscreen = false,
                  sdev_prop         = 0.1,
                  sdev_pos          = [10.0, 10.0],
                  pnorm             = 2,
                  debug             = false,
                  fdataname         = "2Dtest_smooth"
                  )
    # this is where we might store the centroids during MARE2DEM inversions
    opt_EM_in  = MCMC_Driver.EMoptions(sd=sd)

    # calculate the RMS on the data (how much noise have we added to the sampled image points)
    ImageRegression.calc_simple_RMS(d, img, opt_in, opt_EM_in, sd)
end

# let's dive into the main MCMC loop
if rank == 0
    @time MCMC_Driver.main(opt_in, d, Tmax, nsamples, opt_EM_in, myMPIparams)
elseif in(rank,WorkerCaptains) && rank != 0
    MCMC_Driver.main(opt_in, d, Tmax, nsamples, opt_EM_in, myMPIparams)
end


# let's wrap things up
MPI.Finalize()
