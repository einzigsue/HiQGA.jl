
# @mpi_do manager begin
    # do stuff in here
#end

function load_m2d(filenameroot,nprocsperchain)

    comm  = MPI.COMM_WORLD
    rank  = MPI.Comm_rank(comm)
    nproc = MPI.Comm_size(comm)

    team         = Integer(floor(rank/nprocsperchain))
    comm_team    = MPI.Comm_split(comm,team,rank)
    rank_team    = MPI.Comm_rank(comm_team)
    # # define the worker captains (the "leader" threads of each team)
    # WorkerCaptains = Array{Int64,1}(0:nprocsperchain:nproc-1)

    # filenameroot = "prism1.0" # MARE2DEM model files to test (i.e. without the .resitivity)
    bquiet       = true       # set to true to turn off all MARE2DEM print statements

    if rank_team == 0

        println("Reading in MARE2DEM files for base filename ",filenameroot)

        M2d = Mare2dem.load(filenameroot,bquiet)

        # #print field names:
        # println("field names in returned struct:")
        # for fname in fieldnames(typeof(M2d))
        #     println("$fname ")
        # end
        # println(" ")

    else
        M2d = Mare2dem.load_dummy(filenameroot,bquiet)
    end

    MPI.Barrier(comm)

end


function m2d_fwd(nprocperchain,M2d)

    comm  = MPI.COMM_WORLD
    rank  = MPI.Comm_rank(comm)
    nproc = MPI.Comm_size(comm)

    team         = Integer(floor(rank/nprocperchain))
    comm_team    = MPI.Comm_split(comm,team,rank)
    rank_team    = MPI.Comm_rank(comm_team)

    bquiet       = true       # set to true to turn off all MARE2DEM print statements

    println("from rank $rank: I'm rank $rank_team on team $team: I'm here!")

    # Set input values:
    if rank_team == 0
        log10rhofree =  1.0 + zeros(size(M2d.log10rhofree))
        savefileroot = "test_fwd_output_julia"
        iteration    = 42  # iteration number (typically). Gets appended as <savefileroot>.<iterationls>.resp
        nd           = length(M2d.data) # needed to allocate space for output response array
        println("num of data: $nd; first few of log10rhofree: $(log10rhofree[1:5])")
    else
        # # input variables only needs to be defined on rank 0 process. other procs get dummy values
        log10rhofree = 0.
        savefileroot = " "
        iteration    = 0.
        nd           = 0
    end

    MPI.Barrier(comm)
    println("Rank $rank, computing forward response!")
    response, chi2misfit = Mare2dem.forward(nprocperteam,log10rhofree,savefileroot,iteration,nd,bquiet)

    println("rank $rank, reporting for duty at the barrier!")
    MPI.Barrier(comm)

end
