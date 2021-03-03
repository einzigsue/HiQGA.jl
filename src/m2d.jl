


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

    println("from rank $rank: I'm rank $rank_team on team $team: I'm here!")


    if rank_team == 0

        println("Reading in MARE2DEM files for base filename ",filenameroot)
        println("My rank is $rank, my team is $team")

        M2d = Mare2dem.load(filenameroot,bquiet)

        # print field names:
        println("field names in returned struct:")
        for fname in fieldnames(typeof(M2d))
            println("$fname ")
        end
        println(" ")

    end

    println("rank $rank, reporting for duting at the barrier!")
    MPI.Barrier(comm)

end
