any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
import dummy_module
using MPI, Distributed, dummy_module

mgr = MPI.start_main_loop(MPI.MPI_TRANSPORT_ALL)

a = rand(5,1)
c = do_a_sum(a)

nchains = 2

rmprocs(workers())

addprocs(nchains)
@info "workers are $(workers())"
@everywhere any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
@everywhere using Distributed, dummy_module

work_saved = Array{Future, 1}(undef, nworkers())

@sync for(idx, pid) in enumerate(workers())
    b = rand(5,1)
    work_saved[idx] = @spawnat pid do_a_sum(b)
    c = @spawnat pid do_a_sum(b)
    d = fetch(work_saved[idx])[1]
    println("happy $d")
end

rmprocs(workers())

MPI.stop_main_loop(mgr)
