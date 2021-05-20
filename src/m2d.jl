
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

    return M2d

end


function m2d_fwd(nprocperchain,M2d)

    comm  = MPI.COMM_WORLD
    rank  = MPI.Comm_rank(comm)
    nproc = MPI.Comm_size(comm)

    team         = Integer(floor(rank/nprocperchain))
    comm_team    = MPI.Comm_split(comm,team,rank)
    rank_team    = MPI.Comm_rank(comm_team)

    bquiet       = true       # set to true to turn off all MARE2DEM print statements

    # println("from rank $rank: I'm rank $rank_team on team $team: I'm here!")
    # println("from rank $rank: first of M2d.log10rhofree = $(M2d.log10rhofree[1])")

    # Set input values:
    if rank_team == 0
        # println("from team rank $rank_team: I'm inside the if statement")
        # println("$(M2d.log10rhofree)")
        tmp = size(M2d.log10rhofree)
        # println("tmp = $tmp")
        log10rhofree =  M2d.log10rhofree
        # println("first of log10rhofree = $(log10rhofree[1])")
        savefileroot = "test_fwd_output_julia"
        # println("savefileroot: $savefileroot")
        iteration    = 42  # iteration number (typically). Gets appended as <savefileroot>.<iterationls>.resp
        # println("iteration: $iteration")
        nd           = length(M2d.data) # needed to allocate space for output response array
        # println("num of data: $nd; first few of log10rhofree: $(log10rhofree[1:5])")
    else
        # # input variables only needs to be defined on rank 0 process. other procs get dummy values
        log10rhofree = M2d.log10rhofree
        savefileroot = " "
        iteration    = 0.
        nd           = 0
    end

    MPI.Barrier(comm)
    # println("Rank $rank, computing forward response! (nprocperchain = $nprocperchain)")
    if rank_team == 1
        # println("nd and iteration: $nd and $iteration")
        # println("bquiet and savefileroot: $bquiet and $savefileroot")
        # println("log10rhofree and nprocperchain: $(log10rhofree) and $nprocperchain")
    end
    response, chi2misfit = Mare2dem.forward(nprocperchain,log10rhofree,savefileroot,iteration,nd,bquiet)

    # println("rank $rank, reporting for duty at the barrier!")
    MPI.Barrier(comm)

    return (response,chi2misfit)

end
