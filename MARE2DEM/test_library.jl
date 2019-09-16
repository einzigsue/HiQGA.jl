import MPI

function main()

    # always initiate MPI with this:
    MPI.Init()

    comm = MPI.COMM_WORLD           # the default MPI communicator, which has all processors invoked with mpirun
    rank = MPI.Comm_rank(comm)     # the id (index) of this processor in communicator comm
    nproc = MPI.Comm_size(comm)     # the total number of processors in the communicator group comm


    #
    # Create local communicators for the PT threads
    #
    MPI.Barrier(comm)

    nProcPerTeam = 4
    team = Integer(floor(rank/nProcPerTeam))

    println(" rank, ptThread: $rank, $team ")


    comm_team = MPI.Comm_split(comm,team,rank)

    rank_team = MPI.Comm_rank(comm_team)

    MPI.Barrier(comm)
    println(" global rank, team rank: $rank, $rank_team ")

    # Test 1: load in file Demo.2.resistivity:
    if rank == 0
        nFreeRegions = 0
        nParamsPerRegion = 0
        nFreeRegions_ref = Ref{Int64}(nFreeRegions)
        nParamsPerRegion_ref = Ref{Int64}(nParamsPerRegion)
        inputResistivityFile = "Demo.2"
        #inputResistivityFile = "simple"
        fileLength = length(inputResistivityFile)
        println("file name length (in Julia): $fileLength")
        ccall( (:mare2dem_load_files_Test, "./mare2dem_lib.o"), Int64,
        (Ref{Int64},Ref{Int64}),nFreeRegions_ref,nParamsPerRegion_ref)
        #ccall( (:mare2dem_load_files, "./mare2dem_lib.o"), Int64,
        #(Ref{String},Ref{Int64},Ref{Int64},Ref{Int64}),inputResistivityFile,fileLength,nFreeRegions_ref,nParamsPerRegion_ref)
        nFreeRegions = nFreeRegions_ref.x
        nParamsPerRegion = nParamsPerRegion_ref.x
        println(" nFreeRegions: $nFreeRegions")
        println(" nParamsPerRegion: $nParamsPerRegion")
    end

    #Test 2: get free parameters and parameter centroids
    if rank == 0
        # allocate the arrays we will pass to mare2dem_get_params
        nLog10RhoFree = Array{Float64}(undef, 0)
        freeCentroids = Array{Float64}(undef, 0, 0)
        ccall( (:mare2dem_get_params, "./mare2dem_lib.o"), Int64,
        (Ref{Float64},Ref{Float64}),nLog10RhoFree,freeCentroids)
        println("num of free params and centroids returned")
    end



    # always end MPI with this:
    #println(" before finalize")
    MPI.Finalize()
    #println(" after finalize")


end

function do_hello(comm,rank,nproc)
    println("Hello world, I am $rank of $nproc")
    MPI.Barrier(comm)               # wait here until all processors reach this command (not really needed, I think)
end





main()
