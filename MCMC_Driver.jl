module MCMC_Driver
using TransD_GP, Distributed, DistributedArrays,
     PyPlot, LinearAlgebra, Formatting, MPI, Statistics

mutable struct EMoptions
    sd      :: Float64
    MLnoise :: Bool
end

mutable struct MPIparams
    #comm            :: Int64
    rank            :: Int64
    nprocs          :: Int64
    team            :: Int64
    teamRank        :: Int64
    WCindex         :: Int64
    WorkerCaptains  :: Array{Int64,1}
end

function MPIparams(;
    #comm,
    rank,
    nprocs,
    team,
    teamRank,
    WCindex,
    WorkerCaptains
    )
    MPIparams(rank,nprocs,team,teamRank,WCindex,WorkerCaptains)
end

function EMoptions(;
            sd = 0.0,
            MLnoise = true)
    @assert sd != 0.0
    EMoptions(sd, MLnoise)
end

struct Tpointer
    fp   :: IOStream
    fstr :: String
end

function get_misfit(m::TransD_GP.Model, d::AbstractArray, opt::TransD_GP.Options, opt_EM::EMoptions)
    chi2by2 = 0.0
    select = .!isnan.(d[:])
    if !opt.debug
        r = m.fstar[select] - d[select]
        if !opt_EM.MLnoise
            chi2by2 = 0.5*norm(r/opt_EM.sd)^2
        else
            N = sum(select)
            chi2by2 = 0.5*N*log(norm(r)^2)
        end
    end
    return chi2by2
end

function mh_step!(m::TransD_GP.Model, d::AbstractArray,
    opt::TransD_GP.Options, stat::TransD_GP.Stats,
    Temp::Float64, movetype::Int, current_misfit::Float64, opt_EM::EMoptions)

    new_misfit = get_misfit(m, d, opt, opt_EM)
    logalpha = (current_misfit - new_misfit)/Temp
    if log(rand()) < logalpha
        current_misfit = new_misfit
        stat.accepted_moves[movetype] += 1
    else
        TransD_GP.undo_move!(movetype, m, opt)
    end
end

function do_mcmc_step(m::TransD_GP.Model, opt::TransD_GP.Options, stat::TransD_GP.Stats,
    current_misfit::Float64, d::AbstractArray,
    Temp::Float64, isample::Int, opt_EM::EMoptions, wp::TransD_GP.Writepointers,
    myMPIparams::MCMC_Driver.MPIparams)

    # select move and do it
    movetype, priorviolate = TransD_GP.do_move!(m, opt, stat)

    println("worker $(myMPIparams.rank) did move #$movetype, was it successful? $priorviolate")

    if !priorviolate
        mh_step!(m, d, opt, stat, Temp, movetype, current_misfit, opt_EM)
    end

    # acceptance stats
    TransD_GP.get_acceptance_stats!(isample, opt, stat)

    # write models
    writemodel = false
    abs(Temp-1.0) < 1e-12 && (writemodel = true)
    TransD_GP.write_history(isample, opt, m, current_misfit, stat, wp, writemodel)

    current_misfit = 0.0
    return current_misfit
end

#=function do_mcmc_step(m::DArray{TransD_GP.Model}, opt::DArray{TransD_GP.Options},
    stat::DArray{TransD_GP.Stats}, current_misfit::DArray{Array{Float64, 1}},
    d::AbstractArray, Temp::Float64, isample::Int, opt_EM::DArray{EMoptions},
    wp::DArray{TransD_GP.Writepointers})

    misfit = do_mcmc_step(localpart(m)[1], localpart(opt)[1], localpart(stat)[1],
                            localpart(current_misfit)[1], localpart(d),
                            Temp, isample, localpart(opt_EM)[1], localpart(wp)[1])

end=#

function close_history(wp::DArray)
    @sync for (idx, pid) in enumerate(procs(wp))
        @spawnat pid TransD_GP.close_history(wp[idx])
    end
end

function open_temperature_file(opt_in::TransD_GP.Options, T::Array{Float64, 1})
    fdataname = opt_in.costs_filename[9:end-4]
    fp_temps  = open(fdataname*"_temps.txt", opt_in.history_mode)
    fmt = "{:d} "
    for i = 1:length(T)-1
        fmt = fmt*"{:f} "
    end
    fmt = fmt*"{:f}"
    tpointer = Tpointer(fp_temps, fmt)
end

function write_temperatures(iter::Int, T::Array{Float64, 1}, tpointer::Tpointer, opt_in::TransD_GP.Options)
    if (mod(iter-1, opt_in.save_freq) == 0 || iter == 1)
        printfmtln(tpointer.fp, tpointer.fstr, iter, T...)
        flush(tpointer.fp)
    end
end

function close_temperature_file(fp::IOStream)
    close(fp)
end

function init_chain_darrays(opt_in::TransD_GP.Options, opt_EM_in::EMoptions, d_in::AbstractArray,
    myMPIparams::MCMC_Driver.MPIparams)

    costs_filename = "misfits_"*opt_in.fdataname
    fstar_filename = "models_"*opt_in.fdataname
    x_ftrain_filename = "points_"*opt_in.fdataname

    idx = myMPIparams.team + 1

    opt_in.costs_filename    = "misfits_"*opt_in.fdataname*"_$idx.bin"
    opt_in.fstar_filename    = "models_"*opt_in.fdataname*"_$idx.bin"
    opt_in.x_ftrain_filename = "points_"*opt_in.fdataname*"_$idx.bin"

    #println("my team: $(myMPIparams.team), $(opt_in.costs_filename)")

    m = TransD_GP.init(opt_in)
    stat = TransD_GP.Stats()
    current_misfit = get_misfit(m, d_in, opt_in, opt_EM_in)
    wp = TransD_GP.open_history(opt_in)

    return m, stat, current_misfit, wp

end

function main(opt_in::TransD_GP.Options, din::AbstractArray, Tmax::Float64, nsamples::Int, opt_EM_in::EMoptions,
    myMPIparams::MCMC_Driver.MPIparams)
    # reestablish these in this new scope, rather than pass them in as arguments
    #=comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    nprocs = MPI.Comm_size(comm)=#

    #println("$(myMPIparams.WorkerCaptains)")

    if myMPIparams.rank == 0
        # establish temperature ladder for PT
        T = 10.0.^range(0, stop = log10(Tmax), length = myMPIparams.nprocs)
        # open the file for storing where the temperatures are at each step in MCMC
        fp_temps = open_temperature_file(opt_in::TransD_GP.Options, T)
        # this array will hold the aggregated misfits across all chains
        misfitAll = zeros(length(T),1)
        # for timing the algorithm
        t2 = time()
    end

    if in(myMPIparams.rank, myMPIparams.WorkerCaptains)
        m, stat, misfit, wp = init_chain_darrays(opt_in, opt_EM_in, din[:], myMPIparams)
        println("my rank is: $(myMPIparams.rank), misfit: $misfit")
        println("$(myMPIparams.WorkerCaptains)")
    end

    comm = MPI.COMM_WORLD
    MPI.Barrier(comm)

    #main MCMC loop
    for isample = 1:nsamples

        #do PT chain swaps (actually, we're just swapping temperatures)
        #first, get the misfits from the worker captains
        if myMPIparams.rank == 0
            for (idx,widx) in enumerate(myMPIparams.WorkerCaptains)
                if widx != 0
                    holder = zeros(1,1)
                    println("preparing to receive from worker $widx...")
                    MPI.Recv!(holder, widx, 0, comm)
                    println("got message from worker $widx, here is misfit: $holder")
                    misfitAll[idx] = holder[1]
                end
            end
        elseif in(myMPIparams.rank,myMPIparams.WorkerCaptains) && myMPIparams.rank != 0
            send_mesg = misfit
            println("worker $(myMPIparams.rank) preparing to send message $misfit to manager...")
            MPI.Send(send_mesg, 0, 0, comm)
            println("message successfully sent from worker $(myMPIparams.rank)")
        end

        if myMPIparams.rank == 0
            println("message passing over, here is the misfit vector")
            println("$misfitAll")
        end

        MPI.Barrier(comm)
        println(" ")

        #next, perform the temperature swaps
        if myMPIparams.rank == 0
            println(" ")
            println("original T: $T")
            println(" ")
            nchains = length(myMPIparams.WorkerCaptains)
            for ichain = nchains:-1:2
                jchain = rand(1:ichain) #now we have two swap candidates (ichain and jchain)
                #do the swap so long as ichain and jchain aren't the same chain!
                if ichain != jchain
                    #log of the swap probability
                    logalpha = (misfitAll[ichain] - misfitAll[jchain]) * (1.0/T[ichain] - 1.0/T[jchain])
                    if log(rand()) < logalpha
                        T[ichain], T[jchain] = T[jchain], T[ichain]
                    end
                end
            end
            println(" ")
            println("final T: $T")
            println(" ")
            #finally, send the new temperatures to their respective chains
            for (idx,widx) in enumerate(myMPIparams.WorkerCaptains)
                if widx != 0
                    send_mesg = T[idx]
                    println("sending message $send_mesg to worker $widx")
                    MPI.Send(send_mesg, widx, 0, comm)
                elseif widx == 0
                    T_local = zeros(1,1)
                    T_local[1] = T[1]
                end
            end
        elseif in(myMPIparams.rank,myMPIparams.WorkerCaptains) && myMPIparams.rank != 0
            sleep(myMPIparams.rank*0.2)
            println("worker $(myMPIparams.rank) waiting to receive message from manager...")
            T_local = zeros(1,1)
            MPI.Recv!(T_local, 0, 0, comm)
            println("worker $(myMPIparams.rank) received message $T_local from manager")
        end
        #
        #end of parallel tempering swaps!
        #

        MPI.Barrier(comm)
        println(" ")

        #do one MCMC step
        if in(myMPIparams.rank,myMPIparams.WorkerCaptains)
            misfit = do_mcmc_step(m, opt_in, stat, misfit, din, T_local[1], isample, opt_EM_in, wp,
                    myMPIparams)
        end

    end # end of loop over nsamples


end # end of main

function nicenup(g::PyPlot.Figure;fsize=14)
    for ax in gcf().axes
        ax.tick_params("both",labelsize=fsize)
        ax.xaxis.label.set_fontsize(fsize)
        ax.yaxis.label.set_fontsize(fsize)
        ax.title.set_fontsize(fsize)
        if typeof(ax.get_legend_handles_labels()[1]) != Array{Any,1}
            ax.legend(loc="best", fontsize=fsize)
        end
    end
    g.tight_layout()
end

end
