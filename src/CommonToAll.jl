module CommonToAll
using PyPlot, StatsBase, Statistics, Distances, LinearAlgebra,
      DelimitedFiles, ..AbstractOperator, NearestNeighbors, Printf, 
      KernelDensitySJ, KernelDensity, Interpolations

import ..Options, ..OptionsStat, ..OptionsNonstat, ..OptionsNuisance,
       ..history, ..GP.κ, ..calcfstar!, ..AbstractOperator.Sounding

export trimxft, assembleTat1, gettargtemps, checkns, getchi2forall, nicenup, plotconv,
        plot_posterior, make1Dhist, make1Dhists, setupz, zcontinue, makezρ, plotdepthtransforms,
        unwrap, getn, geomprogdepth, assemblemodelsatT, getstats, gethimage,
        assemblenuisancesatT, makenuisancehists, stretchexists,
        makegrid, whichislast, makesummarygrid, makearray, plotNEWSlabels, 
        plotprofile, gridpoints, splitsoundingsbyline, dfn2hdr, getgdfprefix, 
        pairinteractionplot, flipline, summaryconductivity, plotsummarygrids1, getVE

# Kernel Density stuff
abstract type KDEtype end
struct SJ <: KDEtype end
struct LSCV <: KDEtype end

function trimxft(opt::Options, burninfrac::Float64, temperaturenum::Int)
    x_ft = assembleTat1(opt, :x_ftrain, burninfrac=burninfrac, temperaturenum=temperaturenum)
    n = assembleTat1(opt, :nodes, burninfrac=burninfrac, temperaturenum=temperaturenum)
    x, ft = zeros(size(opt.xall, 1), sum(n)), zeros(size(opt.fbounds, 1), sum(n))
    nlast = 0
    for (i, xft) in enumerate(x_ft)
        x[:,nlast+1:nlast+n[i]]  = xft[1:size(opt.xall, 1), 1:n[i]]
        ft[:,nlast+1:nlast+n[i]] = xft[size(opt.xall, 1)+1:end, 1:n[i]]
        nlast += n[i]
    end
    x, ft, n
end

function assembleTat1(optin::Options, stat::Symbol; burninfrac=0.5, temperaturenum=-1)
    @assert temperaturenum!=-1 "Please explicitly specify which temperature number"
    if stat == :fstar && temperaturenum!= 1
        return assemblemodelsatT(optin, burninfrac=burninfrac, temperaturenum=temperaturenum)
    end
    isns = checkns(optin)
    @assert 0.0<=burninfrac<=1.0
    Tacrosschains = gettargtemps(optin)
    temps = sort(unique(Tacrosschains))
    niters = size(Tacrosschains, 1)
    start = round(Int, niters*burninfrac)
    start == 0 && (start = 1)
    ttarget = temps[temperaturenum]
    nmodels = sum((Tacrosschains[start:end,:] .== ttarget))
    if stat == :fstar || stat == :x_ftrain
        mat1 = Array{Array{Float64}, 1}(undef, nmodels)
    else
        mat1 = Array{Real, 1}(undef, nmodels)
    end
    opt = deepcopy(optin)
    imodel = 0
    opt.fstar_filename = "models_"*opt.fdataname*isns*".bin"
    opt.x_ftrain_filename = "points_"*opt.fdataname*isns*".bin"
    opt.costs_filename = "misfits_"*opt.fdataname*isns*".bin"
    chain_idx = nothing
    if isfile(opt.fstar_filename)
        chain_idx = 1
    end
    for ichain in 1:size(Tacrosschains, 2)
        if isnothing(chain_idx)
            opt.fstar_filename = "models_"*opt.fdataname*isns*"_$ichain.bin"
            opt.x_ftrain_filename = "points_"*opt.fdataname*isns*"_$ichain.bin"
            opt.costs_filename = "misfits_"*opt.fdataname*isns*"_$ichain.bin"
        else
            chain_idx = ichain
        end
        if stat == :fstar
            at1idx = findall(Tacrosschains[:,ichain].==ttarget) .>= start
        else
            at1idx = Tacrosschains[:,ichain].==ttarget
            at1idx[1:start-1] .= false
        end
        ninchain = sum(at1idx)
        @info "chain $ichain has $ninchain models"
        ninchain == 0 && continue
        mat1[imodel+1:imodel+ninchain] .= history(opt, stat=stat, chain_idx=chain_idx)[at1idx]
        imodel += ninchain
    end
    iters = history(opt, stat=:iter, chain_idx=chain_idx)
    @info "obtained models $(iters[start]) to $(iters[end]) at T=$ttarget"
    mat1
end

function assemblemodelsatT(opt::OptionsStat; burninfrac=0.9, temperaturenum=-1)
    @assert temperaturenum!=-1 "Please explicitly specify which temperature number"
    x_ft = assembleTat1(opt, :x_ftrain, burninfrac=burninfrac, temperaturenum=temperaturenum)
    n = assembleTat1(opt, :nodes, burninfrac=burninfrac, temperaturenum=temperaturenum)
    nmodels = length(n)
    matT = Array{Array{Float64}, 1}(undef, nmodels)
    K_y = zeros(opt.nmax, opt.nmax)
    Kstar = zeros(Float64, size(opt.xall,2), opt.nmax)
    xtest = opt.xall
    for imodel = 1:nmodels
        xtrain = @view x_ft[imodel][1:size(opt.xall, 1), 1:n[imodel]]
        ftrain = x_ft[imodel][size(opt.xall, 1)+1:end, :]
        ky = view(K_y, 1:n[imodel], 1:n[imodel])
        map!(x->κ(opt.K, x),ky,pairwise(WeightedEuclidean(1 ./opt.λ² ), xtrain, dims=2))
        K_y[diagind(K_y)] .+= opt.δ^2
        ks = view(Kstar, :, 1:n[imodel])
        map!(x->κ(opt.K, x),Kstar,pairwise(WeightedEuclidean(1 ./opt.λ² ), xtest, xtrain, dims=2))
        matT[imodel] = zeros(size(opt.fbounds, 1), size(opt.xall, 2))
        fstar = matT[imodel]
        calcfstar!(fstar, ftrain, opt, K_y, Kstar, n[imodel])
    end
    matT
end

function assemblenuisancesatT(optn::OptionsNuisance;
    burninfrac = 0.5, temperaturenum = -1)
    @assert temperaturenum != -1 "Please specify temperature idx explicitly"
    @assert 0.0 <= burninfrac < 1.0
    Tacrosschains = gettargtemps(optn)
    #this is probably insanely inefficient
    #Θ(niters*nchains) to run the unique() call.
    #however, we need to iterate over the entire array
    #at some point to build the ensemble so it's
    #fine. Possible improvement: store temperatures
    #as the "temperature number" to start with and keep
    #a separate mapping from temperature index to the float
    #value. This means our sortedTs is just 1:nchains and
    #also avoids any issues with equality testing of floats.
    sortedTs = sort(unique(Tacrosschains))
    niters = size(Tacrosschains,1)
    #this will never give a bounds error
    #because of the assert above
    firsti = round(Int, niters*burninfrac)
    firsti == 0 && (start = 1)
    # firsti = 1 + floor(Int, niters*burninfrac)
    ttarg = sortedTs[temperaturenum]
    nmodels = sum(Tacrosschains[firsti:end,:] .== ttarg)
    fdataname = optn.fdataname
    #drop iteration number
    mvals = zeros(nmodels,optn.nnu)
    imodel = 0
    for ichain = 1:size(Tacrosschains,2)
        at1idx = Tacrosschains[:,ichain].==ttarg
        at1idx[1:firsti-1] .= false
        ninchain = sum(at1idx)
        @info "chain $ichain has $ninchain models"
        ninchain == 0 && continue
        vals_filename = "values_nuisance_"*fdataname*".bin"
        if isfile(vals_filename)
            nraw = readdlm(vals_filename, ' ', String)
            cids = parse.(Int, nraw[:,1])
            ndat = parse.(Float64, nraw[cids .== ichain, 3:end])
        else
            vals_filename = "values_nuisance_"*fdataname*"$ichain.bin"
            ndat = readdlm(vals_filename, ' ', Float64)[:,2:end]
        end
        mvals[imodel+1:imodel+ninchain,:] .= ndat[at1idx,:]
        imodel += ninchain
    end
    mvals
end

function getnchains(costs_filename)
    if isfile(costs_filename*".bin")
        data = readdlm(costs_filename * ".bin", String)
        chids = parse.(Int, data[:,1])
        return maximum(chids)
    end

    c = 0
    r = Regex(costs_filename)
    for fname in readdir(pwd())
           c += !isa(match(r, fname), Nothing)
    end
    c
end

function gettargtemps(opt_in::Options)
    isns = checkns(opt_in)
    opt = deepcopy(opt_in)
    costs_filename = "misfits_"*opt.fdataname*isns
    nchains = getnchains(costs_filename)
    @info "Number of chains is $nchains"
    # now look at any chain to get how many iterations
    multichainfile = nothing
    if isfile(costs_filename*".bin")
        opt.costs_filename = costs_filename*".bin"
        multichainfile = 1 # set if all chains are in one file
    else
        opt.costs_filename    = costs_filename*"_1.bin"
    end
    iters          = history(opt, stat=:iter, chain_idx=multichainfile)
    niters         = length(iters)
    @info "McMC has run for $(iters[end]) iterations"
    # then create arrays of unsorted by temperature T
    Tacrosschains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    for ichain in 1:nchains
        if isnothing(multichainfile)
            opt.costs_filename = costs_filename*"_$ichain.bin"
            Tacrosschains[:,ichain] = history(opt, stat=:T)
        else
            Tacrosschains[:,ichain] = history(opt, stat=:T, chain_idx=ichain)
        end
    end
    Tacrosschains
end

function gettargtemps(optn_in::OptionsNuisance)
    optn = deepcopy(optn_in)
    costs_filename = "misfits_"*optn.fdataname*"nuisance"
    nchains = getnchains(costs_filename)
    @info "Number of chains is $nchains"
    multichainfile = nothing
    if isfile(costs_filename*".bin")
        optn.costs_filename = costs_filename*".bin"
        multichainfile = 1 # set if all chains are in one file
    else
        optn.costs_filename    = costs_filename*"_1.bin"
    end
    iters = history(optn, stat=:iter, chain_idx=multichainfile)
    niters = length(iters)
    @info "MCMC has run for $(iters[end]) iterations"
    Tacrosschains = zeros(Float64, niters, nchains)
    for ichain in 1:nchains
        if isnothing(multichainfile)
            optn.costs_filename = costs_filename*"_$ichain.bin"
            Tacrosschains[:,ichain] = history(optn, stat=:T)
        else
            Tacrosschains[:,ichain] = history(optn, stat=:T, chain_idx=ichain)
        end
    end
    Tacrosschains
end

function getstats(optin::Options;
                  figsize=(5,6), fontsize=12,
                  nxticks=5, nyticks=5, alpha=0.6)
    isns = checkns(optin)
    opt = deepcopy(optin)
    costs_filename = "misfits_"*opt.fdataname*isns
    nchains = getnchains(costs_filename)
    chains = 1:nchains
    @info "Number of chains is $nchains"
    multichainfile = nothing
    if isfile(costs_filename*".bin")
        opt.costs_filename = costs_filename*".bin"
        multichainfile = 1
    else
        opt.costs_filename    = costs_filename*"_1.bin"
    end
    iters          = history(opt, stat=:iter, chain_idx=multichainfile)
    statnames = [:acceptanceRateBirth, :acceptanceRateDeath,
                 :acceptanceRatePosition, :acceptanceRateProperty, :acceptanceRateDC]
    f,ax = plt.subplots(length(statnames), 1,
                        sharex=true, sharey=true, figsize=figsize)
    maxar = 0
    for (ichain, chain) in enumerate(chains)
        if isnothing(multichainfile)
            opt.costs_filename = costs_filename*"_$chain.bin"
            chain_idx = nothing
        else
            chain_idx = ichain
        end
        for (istat, statname) in enumerate(statnames)
            ar = history(opt, stat=statname, chain_idx=chain_idx)
            mx = maximum(ar[.!isnan.(ar)])
            mx > maxar && (maxar = mx)
            ax[istat].plot(iters, ar, alpha=alpha)
            if ichain == nchains
                ax[istat].set_title("$statname "*isns)
                ax[istat].set_ylabel("acc. %")
                if istat == length(statnames)
                    ax[istat].set_xlabel("mcmc step no.")
                    ax[istat].set_xticks(LinRange(iters[1], iters[end], nxticks))
                    ax[istat].set_xlim(extrema(iters))
                    ax[istat].set_yticks(LinRange(0, maxar, nyticks))
                    ax[istat].set_ylim(0, maxar)
                end
            end
        end
    end
    nicenup(f, fsize=fontsize)
end

function getstats(optin::OptionsNuisance;
                  figsize=(5,6), fontsize=12,
                  nxticks=5, nyticks=5, alpha=0.6)
    opt = deepcopy(optin)
    costs_filename = "misfits_"*opt.fdataname*"nuisance"
    nchains = getnchains(costs_filename)
    chains = 1:nchains
    @info "Number of chains is $nchains"
    multichainfile = nothing
    if isfile(costs_filename*".bin")
        opt.costs_filename = costs_filename*".bin"
        multichainfile = 1
    else
        opt.costs_filename = costs_filename*"_1.bin"
    end
    iters          = history(opt, stat=:iter, chain_idx=multichainfile)
    statname = :acceptanceRate
    f,ax = plt.subplots(length(optin.idxnotzero), 1, sharex=true, sharey=true, figsize=figsize)
    maxar = 0
    for (ichain, chain) in enumerate(chains)
        if isnothing(multichainfile)
            opt.costs_filename = costs_filename*"_$chain.bin"
            chain_idx = nothing
        else
            chain_idx = ichain
        end
        ar = history(opt, stat=statname, chain_idx=chain_idx)[:,optin.idxnotzero]
        for (i, idx) in enumerate(optin.idxnotzero)
            mx = maximum(ar[.!isnan.(ar)])
            mx > maxar && (maxar = mx)
            ax[i].plot(iters, ar[:,i], alpha=alpha)
            if ichain == nchains
                ax[i].set_title("Nuisance #$idx")
                ax[i].set_ylabel("acc. %")
                if i == length(optin.idxnotzero)
                    ax[i].set_xlabel("mcmc step no.")
                    ax[i].set_xticks(LinRange(iters[1], iters[end], nxticks))
                    ax[i].set_xlim(extrema(iters))
                end
            end
        end
    end
    nicenup(f, fsize=fontsize)
end

function checkns(optin::Options)
    isns = typeof(optin) == OptionsNonstat
    @info "ns is $isns"
    ns = "ns"
    isns || (ns="s")
    return ns
end

function getchi2forall(optn_in::OptionsNuisance;
                       nchains = 1,
                       figsize = (6,4),
                       fsize = 8,
                       alpha = 0.25,
                       nxticks = 0,
                       gridon = false)
    optn = deepcopy(optn_in)
    fdataname = optn.fdataname
    costs_filename = "misfits_"*fdataname*"nuisance"
    if nchains == 1
        nchains = getnchains(costs_filename)
    end
    optn.costs_filename = costs_filename*"_1.bin"
    multichainfile = false
    if isfile(optn.costs_filename)
        iters = history(optn, stat=:iter)
    else
        multichainfile = true
        optn.costs_filename = costs_filename*".bin"
        iters = history(optn, stat=:iter, chain_idx=1)
    end
    niters = length(iters)
    Tacrosschains = zeros(Float64, niters, nchains)
    χ2acrosschains = zeros(Float64, niters, nchains)
    chain_idx = nothing
    for chain in 1:nchains
        if multichainfile
            chain_idx = ichain
        else
            optn.costs_filename = costs_filename*"_$(chain).bin"
        end
        Tacrosschains[:,chain] = history(optn, stat=:T, chain_idx=chain_idx)
        χ2acrosschains[:,chain] = history(optn, stat=:misfit, chain_idx=chain_idx)
    end
    Torder = sort([(i,j) for i=1:niters, j=1:nchains],
                    by = ix->Tacrosschains[ix...],
                    dims = 2)
    χ2sorted = [χ2acrosschains[Torder[i,j]...] for i=1:niters, j=1:nchains]
    fig, ax = subplots(3,1, sharex=true, figsize=figsize)
    ax[1].plot(iters, χ2acrosschains, alpha=alpha)
    ax[1].set_xlim(extrema(iters)...)
    ax[1].set_ylim(0,100)
    ax[1].set_title("unsorted χ^2 misfit")
    ax[2].plot(iters, Tacrosschains, alpha=alpha)
    ax[2].set_title("temperature")
    ax[3].plot(iters, χ2sorted, alpha=alpha)
    ax[3].set_ylim(0,100)
    ax[3].set_title("sorted χ^2 misfit")
    ax[3].set_xlabel("iterations")

    nxticks == 0 || ax[3].set_xticks(iters[1]:div(iters[end],nxticks):iters[end])
    nicenup(fig, fsize = fsize)
end

function getchi2forall(opt_in::Options;
                        nchains         = 1,
                        figsize         = (6,4),
                        fsize           = 8,
                        alpha           = 0.5,
                        nxticks         = 0,
                        gridon          = false,
                        omittemp        = false,
                        hidetitle       = true,
                      )
    # now look at any chain to get how many iterations
    isns = checkns(opt_in)
    opt = deepcopy(opt_in)
    costs_filename = "misfits_"*opt.fdataname*isns
    if nchains == 1 # then actually find out how many chains there are saved
        nchains = getnchains(costs_filename)
    end
    opt.costs_filename    = costs_filename*"_1.bin"
    multichainfile = false
    if isfile(opt.costs_filename)
        iters          = history(opt, stat=:iter)
    else
        opt.costs_filename = costs_filename*".bin"
        iters = history(opt, stat=:iter, chain_idx=1)
        multichainfile=true
    end
    niters         = length(iters)
    # then create arrays of unsorted by temperature T, k, and chi2
    Tacrosschains  = zeros(Float64, niters, nchains)
    kacrosschains  = zeros(Int, niters, nchains)
    X2by2inchains  = zeros(Float64, niters, nchains)
    # get the values into the arrays
    chain_idx = nothing
    for ichain in 1:nchains
        if multichainfile
            chain_idx = ichain
        else
            opt.costs_filename = costs_filename*"_$ichain.bin"
        end
        Tacrosschains[:,ichain] = history(opt, stat=:T, chain_idx=chain_idx)
        kacrosschains[:,ichain] = history(opt, stat=:nodes, chain_idx=chain_idx)
        X2by2inchains[:,ichain] = history(opt, stat=:U, chain_idx=chain_idx)
    end

    if !hidetitle # then we are usually not interested in the temperature sorting of chains    
        f, ax = plt.subplots(3,1, sharex=true, figsize=figsize)
        ax[1].plot(iters, kacrosschains, alpha=alpha)
        ax[1].set_xlim(extrema(iters)...)
        ax[1].set_title(isns*" unsorted by temperature")
        gridon && ax[1].grid()
        ax[1].set_ylabel("# nodes")
        ax[2].plot(iters, X2by2inchains, alpha=alpha)
        gridon && ax[2].grid()
        ax[2].set_ylabel("-Log L")
        gridon && ax[3].grid()
        ax[3].plot(iters, Tacrosschains, alpha=alpha)
        ax[3].set_ylabel("Temperature")
        ax[3].set_xlabel("iterations")
        nxticks == 0 || ax[3].set_xticks(iters[1]:div(iters[end],nxticks):iters[end])
        nicenup(f, fsize=fsize)
    end

    Tunsorted = copy(Tacrosschains)
    for jstep = 1:niters
        sortidx = sortperm(vec(Tacrosschains[jstep,:]))
        X2by2inchains[jstep,:] = X2by2inchains[jstep,sortidx]
        kacrosschains[jstep,:] = kacrosschains[jstep,sortidx]
        Tacrosschains[jstep,:] = Tacrosschains[jstep,sortidx]
    end
    nrows = omittemp ? 2 : 3
    f, ax = plt.subplots(nrows, 1, sharex=true, figsize=figsize)
    nchainsatone = sum(Tacrosschains[1,:] .== 1)
    ax[1].plot(iters, kacrosschains, alpha=alpha)
    ax[1].set_xlim(extrema(iters)...)
    ax[1].set_ylabel("# nuclei")
    !hidetitle && ax[1].set_title(isns*" sorted by temperature")
    ax[1].plot(iters, kacrosschains[:,1:nchainsatone], "k", alpha=alpha)
    gridon && ax[1].grid()
    ax[2].plot(iters, X2by2inchains, alpha=alpha)
    ax[2].plot(iters, X2by2inchains[:,1:nchainsatone], "k", alpha=alpha)
    ax[2].set_ylabel("-Log L")
    gridon && ax[2].grid()
    if !omittemp
        ax[3].plot(iters, Tunsorted, alpha=alpha, color="gray")
        ax[3].set_ylabel("Temperature")
        gridon && ax[3].grid()
    end
    nxticks == 0 || ax[3].set_xticks(iters[1]:div(iters[end],nxticks):iters[end])
    ax[nrows].set_xlabel("iterations")
    nicenup(f, fsize=fsize)

end

function plotconv(optin::Options; burninfrac= 0.5, first=0.1, last=0.5, till=1.0, figsize=(6,2), fontsize=10, nbins=50, doall=false)
    opt = deepcopy(optin)
    iter, k, chi2by2 = map(x->assembleTat1(opt, x; temperaturenum=1, burninfrac), (:iter, :nodes, :U))
    idx = sortperm(iter)
    till = round(Int, idx[end]*till)
    iter, k, chi2by2 = map(x->x[idx][1:till], (iter, k, chi2by2))
    f, ax = plt.subplots(1, 2; figsize)
    # first group, second group
    getconv(ax, iter, chi2by2, k, nbins, opt.nmin, opt.nmax, first, last)
    ax[1].set_xlabel("-ve log likelihood")
    ax[1].set_ylabel("pdf")
    ax[2].set_xlabel("# nuclei")
    ax[2].set_ylabel("pdf")
    nicenup(f, fsize=fontsize)
    nothing
end 

function getconv(ax, iter, chi2by2, k, nbins, nmin, nmax, first, last)
    first, last = map(x->round(Int, length(iter)*x), (first, last))
    plotconv(ax, chi2by2[1:first], k[1:first], nbins, nmin, nmax)
    plotconv(ax, chi2by2[last:end], k[last:end], nbins, nmin, nmax)
end    

function plotconv(ax, chi2by2, k, nbins, nmin, nmax)
    kdefunc_nll = kde_(SJ(), chi2by2)
    nllrange = range(minimum(chi2by2), maximum(chi2by2), nbins)
    krange = nmin-0.5:1:nmax+0.5
    ax[1].plot(nllrange, kdefunc_nll(nllrange))
    ax[2].hist(k, krange, density=true, histtype="step")
end

function geomprogdepth(n, dy, c)
    dy*(1.0-c^n)/(1-c)
end

function getn(z, dy, c)
    log(1 - z/dy*(1-c))/log(c)
end

function unwrap(v, inplace=false)
  # currently assuming an array
  unwrapped = inplace ? v : copy(v)
  for i in 2:length(v)
    while unwrapped[i] - unwrapped[i-1] >= pi
      unwrapped[i] -= 2pi
    end
    while unwrapped[i] - unwrapped[i-1] <= -pi
      unwrapped[i] += 2pi
    end
  end
  return unwrapped
end

function setupz(zstart, extendfrac;n=100, dz=10, showplot=false, forextension=false)
    @assert extendfrac>1
    znrange        = 1.0:n
    zboundaries    = [zstart, zstart .+ geomprogdepth.(znrange, dz, extendfrac)...]
    thickness      = diff(zboundaries)
    zall           = zboundaries[1:end-1] + thickness/2
    znall          = getn.(zall .- zstart, dz, extendfrac)
    showplot && plotdepthtransforms(zall, znall, zboundaries)
    zbreturn = forextension ? zboundaries : zboundaries[1:end-1]
    return zall, znall, zbreturn
end

function zcontinue(;zall=nothing, znall=nothing, zboundaries=nothing, n=nothing,
                extendfrac=nothing, dz=nothing, showplot=false)
    isa(dz, Nothing) && (dz=diff(zboundaries)[end])
    zall_, znall_, zboundaries_ = setupz(zboundaries[end], extendfrac, n=n, dz=dz, showplot=false)
    zall, znall, zboundaries = [zall;zall_], [znall; znall[end] .+ znall_], [zboundaries[1:end-1]; zboundaries_]
    zlast = zboundaries[end] + diff(zboundaries)[end]
    showplot && plotdepthtransforms(zall, znall, [zboundaries;zlast])
    @assert all(diff(zall).>0)
    @assert all(diff(znall).>0)
    @assert all(diff(zboundaries).>0)
    zall, znall, zboundaries
end

function setupz(zstart;n=100, dz=10)
    zboundaries    = range(zstart, step=dz, length=n)
    zall           = zboundaries[1:end-1] .+ dz/2
    znall          = 1.5:1:n-0.5
    plotdepthtransforms(zall, znall, zboundaries)
    return zall, znall, zboundaries[1:end-1]
end

function plotdepthtransforms(zall, znall, zboundaries; fontsize=12)
    thickness = diff(zboundaries)
    f, ax = plt.subplots(2, 2, sharex="col", sharey="row", figsize=(8,5))
    ax[1].stem(zboundaries[1:end-1], zboundaries[1:end-1], markerfmt="")
    ax[1].stem(zall, zall, "k--", markerfmt=" ")
    ax[1].set_ylabel("depth m")
    ax[2].stem(zall,thickness, "k--", markerfmt=" ")
    ax[2].set_ylabel("thickness m")
    ax[2].set_xlabel("depth m")
    ax[2].yaxis.grid(which="major")
    ax[3].plot(znall, zall)
    ax[3].grid()
    ax[4].stem(znall,thickness, "k--", markerfmt=" ")
    ax[4].set_xlabel("depth index")
    ax[4].yaxis.grid(which="major")
    plt.suptitle("Forward model discretization"; fontsize)
    nicenup(f, fsize=fontsize)
end

function makezρ(zboundaries::Array{Float64, 1};
          zfixed      = [-1e6,    0.],
          ρfixed      = [1e12     0.3])
    @assert length(zfixed) == length(ρfixed)
    z = [zfixed..., zboundaries...]
    @assert all(diff(z).>0)
    ρ = zeros(length(z))
    nfixed = length(ρfixed)
    ρ[1:nfixed] .= ρfixed
    z, ρ, nfixed
end

function nicenup(g::PyPlot.Figure;fsize=16, h_pad=nothing)
    for ax in gcf().axes
        ax.tick_params("both",labelsize=fsize)
        ax.xaxis.label.set_fontsize(fsize)
        ax.yaxis.label.set_fontsize(fsize)
        any(keys(ax) .== :zaxis) && ax.zaxis.label.set_fontsize(fsize)
        ax.title.set_fontsize(fsize)

        if typeof(ax.get_legend_handles_labels()[1]) != Array{Any,1}
            ax.legend(loc="best", fontsize=fsize)
        end
    end
    if isnothing(h_pad)
        g.tight_layout()
    else
        g.tight_layout(;h_pad)
    end        
end

function plot_posterior(F::Operator1D,
                        optns::OptionsNonstat,
                        opts::OptionsStat;
                        temperaturenum = 1,
                        nbins = 50,
                        burninfrac=0.5,
                        qp1=0.05,
                        qp2=0.95,
                        vmaxpc=1.0,
                        cmappdf = "inferno",
                        figsize=(12,5),
                        pdfnormalize=false,
                        fsize=14,
                        showlscale1vd=false,
                        istothepow=false,
                        usekde = false,
                        kdetype = SJ(),
                        doplot = true,
                        a = 0.25,
                        lw = 1)
    @assert 0<vmaxpc<=1
    M = assembleTat1(optns, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage_ns, edges_ns, CI_ns, meanimage_ns, 
    meandiffimage_ns, sdslope_ns, = gethimage(F, M, optns; temperaturenum=temperaturenum,
                            nbins=nbins, qp1=qp1, qp2=qp2, istothepow=istothepow, usekde, kdetype,
                            islscale=false, pdfnormalize=pdfnormalize)
    M = assembleTat1(opts, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage, edges, CI, meanimage, 
    meandiffimage, sdslope, = gethimage(F, M, opts; temperaturenum=temperaturenum,
                            nbins=nbins, qp1=qp1, qp2=qp2, istothepow=false, usekde, kdetype,
                            islscale=true, pdfnormalize=pdfnormalize)
    f,ax = plt.subplots(1, 3+convert(Int, showlscale1vd), sharey=true, figsize=figsize)
    xall = opts.xall
    diffs = diff(xall[:])
    xmesh = vcat(xall[1:end-1] - diffs/2, xall[end]-diffs[end]/2, xall[end])
    vmin, vmax = extrema(himage_ns)
    vmax = vmin+vmaxpc*(vmax-vmin)
    im1 = ax[1].pcolormesh(edges_ns[:], xmesh, himage_ns, cmap=cmappdf, vmax=vmax)
    cb1 = colorbar(im1, ax=ax[1])
    cb1.ax.set_xlabel("pdf \nns")
    @show propmin, propmax = getbounds(CI_ns, optns.fbounds)
    ax[1].set_xlim(propmin, propmax)
    ax[2].fill_betweenx(xall[:],meandiffimage_ns[:]-sdslope_ns[:],meandiffimage_ns[:]+sdslope_ns[:],alpha=a, color="gray")
    ax[2].plot(meandiffimage_ns, xall[:], linewidth=lw, color="k")
    ax[2].set_xlabel("slope")

    vmin, vmax = extrema(himage)
    vmax = vmin+vmaxpc*(vmax-vmin)
    im3 = ax[3].pcolormesh(edges[:], xmesh, himage, cmap=cmappdf, vmax=vmax)
    ax[3].set_ylim(extrema(xall)...)
    propmin, propmax = getbounds(CI, opts.fbounds)
    ax[3].set_xlim(propmin, propmax)
    ax[3].invert_yaxis()
    cb3 = colorbar(im3, ax=ax[3])
    cb3.ax.set_xlabel("pdf \nstationary")

    if showlscale1vd
        ax[4].plot(meandiffimage[:], xall[:], linewidth=2, color="g", alpha=0.5)
        ax[4].fill_betweenx(xall[:], meandiffimage[:]-sdslope[:],meandiffimage[:]+sdslope[:], alpha=.5)
        ax[4].set_xlabel("slope")
    end

    ax[1].plot(CI_ns, xall[:], linewidth=lw, color="w", alpha=a)
    # ax[1].plot(CI_ns, xall[:], linewidth=1, color="c", linestyle="--", alpha=0.5)
    # ax[1].plot(meanimage_ns[:], xall[:], linewidth=2, color="w", alpha=0.5)
    # ax[1].plot(meanimage_ns[:], xall[:], linewidth=2, color="k", linestyle="--", alpha=0.5)
    ax[3].plot(CI, xall[:], linewidth=lw, color="w", alpha=a)
    # ax[3].plot(CI, xall[:], linewidth=1, color="c", linestyle="--", alpha=0.5)
    # ax[3].plot(meanimage[:], xall[:], linewidth=2, color="w", alpha=0.5)
    # ax[3].plot(meanimage[:], xall[:], linewidth=2, color="k", linestyle="--", alpha=0.5)
    ax[1].set_xlabel(L"\log_{10} \rho")
    ax[1].set_ylabel("depth (m)")
    ax[3].set_xlabel(L"\log_{10} λ")
    nicenup(f, fsize=fsize)
end

function plot_posterior(F::Operator1D,
                    opt::OptionsStat;
                    useoptfbounds = true,
                    temperaturenum = 1,
                    nbins = 50,
                    burninfrac=0.5,
                    qp1=0.05,
                    qp2=0.95,
                    vmaxpc=1.0,
                    cmappdf = "inferno",
                    figsize=(5,5),
                    pdfnormalize=false,
                    istothepow=false,
                    fsize=14,
                    plotCI = true,
                    plotmean = true,
                    alpha = 1.0,
                    CIcolor = ["w", "k"],
                    meancolor = ["m", "r"],
                    lwidth = 2,
                    pdfclim = nothing,
                    showslope = false,
                    kdetype = SJ(),
                    usekde = false,
                    doplot = true)
    @assert 0<vmaxpc<=1
    
    if useoptfbounds
        rhomin, rhomax = extrema(opt.fbounds)
    else
        rhomin, rhomax = Inf, -Inf
    end    
    M = assembleTat1(opt, :fstar, burninfrac=burninfrac, temperaturenum=temperaturenum)
    himage, edges, CI, meanimage, meandiffimage, sdslope, = gethimage(F, M, opt; temperaturenum,
                nbins, qp1, qp2, istothepow, rhomin, rhomax, usekde, kdetype,
                islscale=false, pdfnormalize=pdfnormalize)

    if doplot
        if showslope
            f, ax = plt.subplots(1,2, sharey=true, figsize=figsize)
        else
            f, ax = plt.subplots(1,1, sharey=true, figsize=figsize, squeeze=false)
        end    
        xall = opt.xall
        diffs = diff(xall[:])
        xmesh = vcat(xall[1:end-1] - diffs/2, xall[end]-diffs[end]/2, xall[end])
        vmin, vmax = extrema(himage)
        vmax = vmin+vmaxpc*(vmax-vmin)
        im1 = ax[1].pcolormesh(edges[:], xmesh, himage, cmap=cmappdf, vmax=vmax)
        ax[1].set_ylim(extrema(xall)...)
        ax[1].invert_yaxis()
        plotCI && ax[1].plot(CI, xall[:], linewidth=lwidth, color=CIcolor[1], alpha=alpha)
        plotCI && ax[1].plot(CI, xall[:], linewidth=lwidth, color=CIcolor[2], linestyle="--", alpha=alpha)
        plotmean && ax[1].plot(meanimage[:], xall[:], linewidth=lwidth, color=meancolor[1], alpha=alpha)
        plotmean && ax[1].plot(meanimage[:], xall[:], linewidth=lwidth, color=meancolor[2], linestyle="--", alpha=alpha)
        ax[1].set_xlabel(L"\log_{10} \rho")
        ax[1].set_ylabel("depth (m)")
        bounds = copy(opt.fbounds)
        if stretchexists(F)
            bounds = [minimum(F.low) maximum(F.low + F.Δ)]
        end
        istothepow && (bounds .= 10 .^ bounds)
        propmin, propmax = getbounds(CI, bounds)
        ax[1].set_xlim(propmin, propmax)
        cb1 = colorbar(im1, ax=ax[1])
        cb1.ax.set_title("pdf")
        if showslope
            ax[2].plot(meandiffimage[:], xall[:], linewidth=2, color="k", linestyle="-")
            zeroside = meandiffimage[:]-sdslope[:]
            zeroside[zeroside .< 0] .= 0
            ax[2].fill_betweenx(xall[:],zeroside,meandiffimage[:]+sdslope[:],alpha=.25)
            ax[2].set_xlabel("mean slope")
        end
        !isnothing(pdfclim) && ax[1].collections[1].set_clim(pdfclim)
        nicenup(f, fsize=fsize)
    end
    CI[:,1], CI[:,2], CI[:,3], meanimage, meandiffimage, sdslope
end

function getbounds(CI, bounds)
    propmin = min(minimum(CI), minimum(bounds))
    propmax = max(maximum(CI), maximum(bounds))
    propmin, propmax
end    

function plot_posterior(operator::Operator1D,
                        optn::OptionsNuisance;
                        temperaturenum = 1,
                        nbins = 50,
                        burninfrac=0.5,
                        figsize=(5,8),
                        pdfnormalize=false,
                        fsize=10,
                        qp1 = 0.05,
                        qp2 = 0.95,
                        labels= nothing,
                        doplot=true)
    hists, CI = makenuisancehists(optn, qp1, qp2, burninfrac = burninfrac,
                                     nbins = nbins, temperaturenum = temperaturenum,
                                 )
    if doplot
        fig,ax = subplots(length(hists), 1, figsize=figsize, squeeze=false)
        for (i, h) = enumerate(hists)
            bwidth, bx, denom = getbinsfromhist(h, pdfnormalize=pdfnormalize)
            ax[i].bar(bx, h.weights./denom, width=bwidth, edgecolor="black")
            !isnothing(labels) && ax[i].set_xlabel(labels[i])
            ax[i].set_ylabel("pdf")
        end
        nicenup(fig, fsize=fsize)
    end
    hists, CI
end

function getbinsfromhist(h, ;pdfnormalize=false)
    bwidth = diff(h.edges[1])
    bx = h.edges[1][1:end-1] + bwidth/2
    denom = pdfnormalize ? maximum(h.weights) : sum(h.weights)*bwidth
    bwidth, bx, denom
end        

function makenuisancehists(optn::OptionsNuisance, qp1, qp2; burninfrac = 0.5, nbins = 50,
    temperaturenum = -1)
    @assert temperaturenum != -1 "Please explicitly specify the temperature index"
    nuisanceatT = assemblenuisancesatT(optn, burninfrac = burninfrac,
                                    temperaturenum = temperaturenum)
    nnu = length(optn.idxnotzero)
    CI = zeros(Float64, nnu, 3)
    h = Vector{Histogram}(undef, nnu)
    for (i, islice) in enumerate(optn.idxnotzero)
        numin, numax = extrema(optn.bounds[islice,:])
        edges = LinRange(numin, numax, nbins+1)
        h[i] = fit(Histogram, nuisanceatT[:,islice], edges)
        CI[i,:] = quantile(nuisanceatT[:,islice], [qp1, 0.5, qp2])
    end
    h, CI
end

function firstderiv(x)
    fd = similar(x)
    fd[1]        =  x[2]     - x[1]
    fd[end]      =  x[end]   - x[end-1]
    fd[2:end-1] .= (x[3:end] - x[1:end-2])/2
    abs.(fd)
    # fd
end

function secondderiv(x)
    sd = similar(x)
    sd[1]        = 2x[1] - 5x[2] + 4x[3] - x[4]
    sd[2:end-1] .= (x[3:end] - x[2:end-1] + x[1:end-2])
    sd[end]      = 2x[end] - 5x[end-1] + 4x[end-2] - x[end-3]
    abs.(sd)
end

function dodensityestimate(usekde::Bool, data, K::KDEtype, edges)
    if usekde
        kdefunc = kde_(K, data)
        kdefunc(0.5(edges[2:end]+edges[1:end-1]))
    else
        w = fit(Histogram, data, edges).weights
        w/(sum(w)*diff(edges)[1])
    end    
end

function stretchexists(F::Operator)
    in(:stretch, fieldnames(typeof(F))) && F.stretch
end

kde_(K::SJ, data) = kde_sj(;data)
kde_(K::LSCV, data) = kde_cv(;data)

abstract type KDEstimator end

# Sheather-Jones plugin
struct kde_sj <: KDEstimator
    data
    bw :: Real
end    

kde_sj(;data=zeros(0)) = kde_sj(data, bwsj(data))

function (foo::kde_sj)(xvals)
    density(foo.data, foo.bw, xvals)
end

# LSCV KDE
struct kde_cv <: KDEstimator
    U
end    

kde_cv(;data=zeros(0)) = kde_cv(kde_lscv(data))

function (foo::kde_cv)(xvals)
    pdf(foo.U, xvals)
end

function gethimage(F::Operator, M::AbstractArray, opt::Options;
                burninfrac = 0.5,
                nbins = 50,
                rhomin=Inf,
                rhomax=-Inf,
                qp1=0.05,
                qp2=0.95,
                islscale = false,
                pdfnormalize=false,
                istothepow = false,
                usekde = false,
                kdetype = SJ(),
                temperaturenum=1)
    T = x->x
    if (rhomin == Inf) && (rhomax == -Inf)
        if !stretchexists(F) # if no stretch
            for (i,mm) in enumerate(M)
                rhomin_mm = minimum(mm)
                rhomax_mm = maximum(mm)
                rhomin_mm < rhomin && (rhomin = rhomin_mm)
                rhomax_mm > rhomax && (rhomax = rhomax_mm)
            end
            if (typeof(opt) == OptionsStat && opt.needλ²fromlog) && islscale
                T = x->0.5log10(x)
            end
        else # there is a stretch
            rhomin = minimum(F.low)
            rhomax = maximum(F.low + F.Δ)
        end        
    else
        @assert rhomin < rhomax
    end
    if istothepow && !islscale
        T = x->10. ^x
    end
    rhomin, rhomax = map(x->T(x), (rhomin, rhomax))    
    edges = LinRange(rhomin, rhomax, nbins+1)
    himage = zeros(Float64, length(M[1]), nbins)
    CI = zeros(Float64, length(M[1]), 3)
    meanimage = zeros(length(M[1]))
    for ilayer=1:length(M[1])
        if !stretchexists(F) # if no stretch
            mthislayer = [T(m[ilayer]) for m in M]
            himage[ilayer,:] = dodensityestimate(usekde, mthislayer, kdetype, edges)
            CI[ilayer,:] = [quantile(mthislayer,(qp1, 0.5, qp2))...]
            meanimage[ilayer] = mean(vec(mthislayer))
        else # there is an affine stretch with depth
            mthislayer = [m[ilayer] for m in M]
            expandedthislayer = T.(F.low[ilayer] .+ mthislayer*F.Δ[ilayer])
            himage[ilayer,:] = dodensityestimate(usekde, expandedthislayer, kdetype, edges)
            CI[ilayer,:] = [quantile(expandedthislayer,(qp1, 0.5, qp2))...]
            meanimage[ilayer] = mean(vec(expandedthislayer))
        end        
        pdfnormalize && (himage[ilayer,:] = himage[ilayer,:]/maximum(himage[ilayer,:]))
    end
    if !stretchexists(F)
        Tm = [T.(m) for m in M]
    else
        Tm = [T.(F.low + vec(m).*F.Δ) for m in M]
    end    
    Mslope = mean(firstderiv.(Tm))
    sdevslope = std(firstderiv.(Tm))
    himage, edges, CI, meanimage, Mslope, sdevslope
end

# plotting codes for 2D sections in AEM

function splitsoundingsbyline(soundings::Array{S, 1}) where S<:Sounding
    alllines = [s.linenum for s in soundings]
    linenum  = unique(alllines)
    nlines   = length(linenum)
    linestartidx = zeros(Int, nlines)
    for i = 1:nlines
        linestartidx[i] = findfirst(alllines .== linenum[i])
    end
    linestartidx
end 

function makegrid(vals::AbstractArray, soundings::Array{S, 1}; donn=false,
    dr=10, zall=[NaN], dz=-1) where S<:Sounding
    @assert all(.!isnan.(zall)) 
    @assert dz>0
    X = [s.X for s in soundings]
    Y = [s.Y for s in soundings]
    R = cumulativelinedist(X,Y)
    if donn
        rr, zz = [r for z in zall, r in R], [z for z in zall, r in R]
        topo = [s.Z for s in soundings]
        zz = topo' .- zz # mAHD
        kdtree = KDTree([rr[:]'; zz[:]'])
        gridr = range(R[1], R[end], step=dr)
        gridz = reverse(range(extrema(zz)..., step=dz))
        rr, zz = [r for z in gridz, r in gridr], [z for z in gridz, r in gridr]
        idxs, = nn(kdtree, [rr[:]'; zz[:]'])
        img = zeros(size(rr))
        for i = 1:length(img)
            img[i] = vals[idxs[i]]
        end
        topofine = gridpoints(R, gridr, topo)
    else
        nodes = ([z for z in zall], [r for r in R])
        itp = extrapolate(interpolate(nodes, vals, Gridded(Linear())), Line()) 
        topo = [s.Z for s in soundings]
        gridr = range(R[1], R[end], step=dr)
        topofine = (interpolate((R,), topo, Gridded(Linear())))(gridr)
        # height = topo in mAHD - depth # mAHD
        minahd = minimum(topo) - maximum(zall)
        maxahd = maximum(topo)
        gridz = reverse(range(minahd, maxahd, step=dz))
        img = [itp(topofine[iy] - x,y) for x in gridz, (iy, y) in enumerate(gridr)]
        zz = [z for z in gridz, r in gridr] 
    end    
    img[zz .>topofine'] .= NaN
    img, gridr, gridz, topofine, R
end

function cumulativelinedist(X,Y)
    dx = diff(X)
    dy = diff(Y)
    R = [0.; cumsum(sqrt.(dx.^2 + dy.^2))]
end    

function gridpoints(R, gridr, points::Array{Float64, 1})
    @assert length(R) == length(points)
    kdtree = KDTree(R[:]')
    idxs, = nn(kdtree, gridr[:]')
    pointsgrid = [points[idxs[i]] for i = 1:length(gridr)]
end

function gridpoints(R, gridr, points::Array{Float64, 2})
    nsound = length(R)
    nprops = size(points, 1)
    @assert nsound == size(points, 2)
    pointsgrid = zeros(length(gridr), nprops)
    for i = 1:nprops
        pointsgrid[:,i] = gridpoints(R, gridr, vec(points[i,:]))
    end
    pointsgrid    
end    

function whichislast(soundings::AbstractArray)
    X = [s.X for s in soundings]
    Y = [s.Y for s in soundings]
    ΔX = X[end] - X[1]
    ΔY = Y[end] - Y[1]
    NSline, EWline = false, false
    if abs(ΔX/ΔY) < 0.05
        NSline = true
    elseif abs(ΔY/ΔX) < 0.05
        EWline = true
    end    
    
    Eislast, Nislast = true, true
    ΔX<0 && (Eislast = false)
    ΔY<0 && (Nislast = false)

    Eislast, Nislast, EWline, NSline
end

function makesummarygrid(soundings, pl, pm, ph, ρmean, vdmean, vddev, zall, dz; dr=10)
    # first flip ρ to σ and so pl and ph are interchanged
    phgrid, gridx, gridz, topofine, R = makegrid(-pl, soundings, zall=zall, dz=dz, dr=dr)
    plgrid,                           = makegrid(-ph, soundings, zall=zall, dz=dz, dr=dr)
    pmgrid,                           = makegrid(-pm, soundings, zall=zall, dz=dz, dr=dr)
    σmeangrid,                        = makegrid(-ρmean, soundings, zall=zall, dz=dz, dr=dr)
    ∇zmeangrid,                       = makegrid(vdmean, soundings, zall=zall, dz=dz, dr=dr)
    ∇zsdgrid,                         = makegrid(vddev, soundings, zall=zall, dz=dz, dr=dr)
    Z = [s.Z for s in soundings]
    phgrid, plgrid, pmgrid, σmeangrid, ∇zmeangrid, ∇zsdgrid, gridx, gridz, topofine, R, Z
end

function makexyzrho(soundings, zall)
    linename = "_line_$(soundings[1].linenum)_summary.txt"
    fnames = ["rho_low", "rho_mid", "rho_hi", "rho_avg"].*linename
    dlow, dmid, dhigh, davg = map(x->readdlm(x), fnames)
    ndepths = length(zall)
    for (i, d) in enumerate([dlow, dmid, dhigh, davg])
        @assert ndims(d) == 2
        @assert size(d, 2) == length(soundings)
        @assert size(d, 1) == ndepths
        xyzrho = makearray(soundings, d, zall)
        writedlm(fnames[i][1:end-4]*"_xyzrho.txt", xyzrho)
    end
end

function makearray(soundings, d, zall)
    outarray = zeros(length(d), 4)
    ndepths = length(zall)
        for (i, s) in enumerate(soundings)
            z = s.Z .- zall # convert depths to mAHD
            x, y = s.X, s.Y
            idx = (i-1)*ndepths+1:i*ndepths
            outarray[idx,:] = [[x y].*ones(ndepths) z d[:,i]]
        end
    outarray
end

function plotNEWSlabels(soundings, gridx, gridz, axarray, 
        x0=nothing, y0=nothing, xend=nothing, yend=nothing;
        preferEright=false, preferNright=false, fontsize=12)
    Eislast, Nislast, EWline, NSline = whichislast(soundings)
    beginpos, endpos = "", ""
    if !any(isnothing.([x0,y0,xend,yend]))
        beginpos = @sprintf(" %.2f", x0/1000)*@sprintf(" %.2f", y0/1000)
        endpos = @sprintf(" %.2f", xend/1000)*@sprintf(" %.2f", yend/1000)
    end
    for s in axarray
        minylim = minimum(s.get_ylim())
        inverted = false
        if (preferNright && !Nislast) || (preferEright && !Eislast) 
            inverted = true
        end 
        
        if Eislast & Nislast # old WN
            dir1, dir2 = "SW", "NE"
            if EWline
                dir1, dir2 = "W", "E"
            elseif NSline
                dir1, dir2 = "S", "N"
            end    
        elseif !Eislast & Nislast # old EN
            dir1, dir2 = "SE", "NW"
            if EWline
                dir1, dir2 = "E", "W"
            elseif NSline
                dir1, dir2 = "S", "N"
            end 
        elseif Eislast & !Nislast # old WS
            dir1, dir2 = "NW", "SE"
            if EWline
                dir1, dir2 = "W", "E"
            elseif NSline
                dir1, dir2 = "N", "S"
            end 
        else # old ES
            dir1, dir2 = "NE", "SW"
            if EWline
                dir1, dir2 = "E", "W"
            elseif NSline
                dir1, dir2 = "N", "S"
            end 
        end
        ha = inverted ? "right" : "left"
        s.text(gridx[1], minylim, dir1*beginpos, backgroundcolor="w", ha = ha; fontsize)
        ha = inverted ? "left" : "right"
        s.text(gridx[end], minylim, dir2*endpos, backgroundcolor="w", ha = ha; fontsize)      
    end
    flipline(preferNright, preferEright, Nislast, Eislast, axarray[end])
end

function flipline(preferNright, preferEright, Nislast, Eislast, ax)
    @assert !(preferEright & preferNright) #both true means don't know what you're doing
    (preferNright && !Nislast) && ax.invert_xaxis()
    (preferEright && !Eislast) && ax.invert_xaxis()
end

function plotprofile(ax, idxs, Z, R)
    for idx in idxs
        ax.plot(R[idx]*[1,1], [ax.get_ylim()[1], Z[idx]], "-w")
        ax.plot(R[idx]*[1,1], [ax.get_ylim()[1], Z[idx]], "--k")
    end
end

function plotsummarygrids1(soundings, meangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname; qp1=0.05, qp2=0.95,
                        figsize=(10,10), fontsize=12, cmap="turbo", vmin=-2, vmax=0.5, 
                        topowidth=2, idx=nothing, omitconvergence=false, useML=false, preferEright=false, preferNright=false,
                        saveplot=false, yl=nothing, dpi=300, showplot=true, showmean=false)
    dr = diff(gridx)[1]
    nrows = omitconvergence ? 5 : 6
    height_ratios = omitconvergence ? [1, 1, 1, 1, 0.1] : [0.4, 1, 1, 1, 1, 0.1]
    if !showmean 
        height_ratios = [height_ratios[1:2]..., height_ratios[4:end]...]
        nrows-=1
    end    
    f, s = plt.subplots(nrows, 1, gridspec_kw=Dict("height_ratios" => height_ratios),
                        figsize=figsize)
    f.suptitle(lname*" Δx=$dr m, Fids: $(length(R))", fontsize=fontsize)
    icol = 1
    if !omitconvergence
        s[icol].plot(R, χ²mean)
        s[icol].plot(R, ones(length(R)), "--k")
        s[icol].fill_between(R, vec(χ²mean-χ²sd), vec(χ²mean+χ²sd), alpha=0.5)
        s[icol].set_ylabel(L"ϕ_d")
        titlestring = useML ? "Max likelihood variance adjustment" : "Data misfit"
        s[icol].set_title(titlestring)
        icol += 1
    end

    summaryconductivity(s, icol, f, soundings, meangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, ; qp1, qp2, fontsize, 
        cmap, vmin, vmax, topowidth, idx, omitconvergence, preferEright, preferNright, yl, showmean)
    
    saveplot && savefig(lname*".png", dpi=dpi)
    showplot || close(f)
end

function summaryconductivity(s, icol, f, soundings, meangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, ; qp1=0.05, qp2=0.95,
    fontsize=12, cmap="turbo", vmin=-2, vmax=0.5, topowidth=2, idx=nothing, omitconvergence=false, preferEright=false, preferNright=false,
    yl=nothing,showmean=false)
    icolstart = icol
    s[icol].imshow(plgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    s[icol].plot(gridx, topofine, linewidth=topowidth, "-k")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    s[icol].set_title("Percentile $(round(Int, 100*qp1)) conductivity")
    s[icol].set_ylabel("Height m")
    omitconvergence || s[icol].sharex(s[icol-1])
    icol += 1
    s[icol].imshow(pmgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    s[icol].plot(gridx, topofine, linewidth=topowidth, "-k")
    s[icol].set_title("Percentile 50 conductivity")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    s[icol].set_ylabel("Height m")
    s[icol].sharex(s[icol-1])
    s[icol].sharey(s[icol-1])
    icol += 1
    if showmean
        s[icol].imshow(meangrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                    extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
        s[icol].plot(gridx, topofine, linewidth=topowidth, "-k")
        s[icol].set_title("Mean conductivity")
        s[icol].set_ylabel("Height m")
        idx == nothing || plotprofile(s[icol], idx, Z, R)
        s[icol].sharex(s[icol-1])
        s[icol].sharey(s[icol-1])
        icol +=1
    end    
    imlast = s[icol].imshow(phgrid, cmap=cmap, aspect="auto", vmax=vmax, vmin = vmin,
                extent=[gridx[1], gridx[end], gridz[end], gridz[1]])
    s[icol].plot(gridx, topofine, linewidth=topowidth, "-k")
    s[icol].set_xlabel("Line distance m", labelpad=0)
    s[icol].set_title("Percentile $(round(Int, 100*qp2)) conductivity")
    s[icol].set_ylabel("Height m")
    idx == nothing || plotprofile(s[icol], idx, Z, R)
    s[icol].sharex(s[icol-1])
    s[icol].sharey(s[icol-1])
    s[icol].set_xlim(extrema(gridx))
    # map(x->x.set_xticklabels([]), s[1:end-2])
    map(x->x.tick_params(labelbottom=false), s[1:end-2])
    # map(x->x.grid(), s[1:end-1])
    isa(yl, Nothing) || s[end-1].set_ylim(yl...)
    plotNEWSlabels(soundings, gridx, gridz, s[icolstart:end-1]; preferEright, preferNright, fontsize)
    cb = f.colorbar(imlast, cax=s[end], orientation="horizontal")
    cb.set_label("Log₁₀ S/m", labelpad=0)
    nicenup(f, fsize=fontsize, h_pad=0)
    label = f._suptitle.get_text()
    VE = round(Int, getVE(s[end-1]))
    f.suptitle(label*", VE=$(VE)X", fontsize=0.8*fontsize)
end

function getVE(ax)
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = abs.(diff([ax.get_ylim()...])[1] / diff([ax.get_xlim()...])[1])
    return disp_ratio / data_ratio
end

function selectwithin1Dinterval(M::AbstractVector, z, 
                                zbounds, vbounds,
                                cond = :mean)
    @assert(any(cond .== [:all, :any, :median, :mean]))
    @assert size(vbounds, 1) == size(zbounds, 1)
    @assert length(z) == length(M[1])
    @assert all(vbounds[:,1] .< vbounds[:,2])
    @assert all(zbounds[:,1] .< zbounds[:,2])
    nconditions = size(vbounds, 1)
    idx = zeros(Bool, length(M), nconditions)
    for (i, m) in enumerate(M)
       for j in 1:nconditions
            idxdepth = zbounds[j,1] .< z .< zbounds[j,2]
            if any(cond .== [:median, :mean])
                if vbounds[j,1] < eval(cond)(m[idxdepth]) < vbounds[j,2]
                    idx[i,j] = true
                end
            else     
                if eval(cond)(vbounds[j,1] .< (m[idxdepth]) .< vbounds[j,2])
                    idx[i,j] = true
                end    
            end
       end
    end
    vec(reduce(&, idx, dims=2))
end

function block1Dvalues(M::AbstractVector, z, zbounds, cond = :median)
    @assert(any(cond .== [:median, :mean]))
    @assert length(z) == length(M[1])
    @assert all(zbounds[:,1] .< zbounds[:,2])
    nconditions = size(zbounds, 1)
    Mblock = zeros(length(M), nconditions)
    for (i, m) in enumerate(M)
        for j in 1:nconditions
            idxdepth = zbounds[j,1] .< z .<= zbounds[j,2]
            Mblock[i,j] = eval(cond)(m[idxdepth])
        end
     end
    Mblock
end    

function correlationplot(M::Array{Float64, 2}; figsize=(5,5), nbins=25)
    f, ax = plt.subplots(size(M, 2), size(M, 2), figsize=figsize, sharex=true, sharey=true)
    nrows = size(M, 2)
    for j = 1:nrows
        for i = j:nrows
            h = normalize(fit(Histogram, (M[:,i], M[:,j]), nbins=nbins))
            ax[nrows*(j-1)+i].pcolormesh(h.edges[2], h.edges[1], h.weights)
        end
    end 
    ax[1].set_aspect(1)   
end    

function plot1Dpatches(ax, zlist, xlist; alpha=0.25, fc="red", ec="blue", lw=2)
    @assert length(zlist) == length(xlist)
    for i in 1:size(zlist, 1)
        x0, z0 = xlist[i,1], zlist[i,1]
        delx = xlist[i,2] - x0
        delz = zlist[i,2] - z0
        ax.add_patch(matplotlib.patches.Rectangle((x0,z0),delx,delz, alpha=alpha, fc=fc, ec=fc))
        ax.add_patch(matplotlib.patches.Rectangle((x0,z0),delx,delz, fill=false, ec=ec, lw=lw, linestyle="--"))
    end    
end

function getclosestidx(Xwell, Ywell, soundings::Vector{S}) where S<: Sounding
    XY = hcat([[s.X, s.Y] for s in soundings]...)
    tree = KDTree(XY)
    idx, = nn(tree, [Xwell;Ywell])
    idx
end    

#function to read the *dfn file and extract the column number and column names as a *.txt file 
function dfn2hdr(dfnfile::String)
    dfn = readlines(dfnfile)
    dfn = dfn[(.!contains.(dfn, "RT=PROJ") .& .!contains.(dfn, "RT=TRNS")
              .& .!contains.(dfn, "RT=COMM"))]

    rectype_rgx = r"RT=([^;]*);"
    name_rgx = r"NAME="
    inc_regex = r":([0-9]+)F" #this will only work with floating point fields

    cumulative_columns = 0  #this will set up a cumulative variable 
    
    fname_hdr = getgdfprefix(dfnfile)*".hdr" # the name of HDR file associated with DFN
    io = open(fname_hdr,"w")
    data_first = false
    for row = dfn
        m = match(rectype_rgx, row)
        isnothing(m) && continue
        data_first |= (m[1] == "DATA")

        inc_match = match(inc_regex, row)  #matching the keyword with the string 
        if isnothing(inc_match)
            inc = 1
            idx2 = split(row, name_rgx)
            idx2_r = first.(split.(idx2[2], ":"))
        else
            inc = parse(Int64, inc_match.captures[1]) #the outcome is int64
            idx2 = split(row, rectype_rgx)
            idx2_r = first.(split.(idx2[2], ":"))
        end

        firstcol = cumulative_columns + 1 + data_first
        lastcol = firstcol + inc - 1
        cumulative_columns += inc
        
        if occursin(";END DEFN", idx2_r)  #type is string 
            idx2_r = first.(split.(idx2_r,";"))
        else
            idx2_r = first.(split.(idx2_r,"END"))
        end

        #start writing to a file here
        if (firstcol != lastcol)
            writedlm(io,[firstcol lastcol idx2_r]) 
        else
            writedlm(io, [firstcol idx2_r])
        end
    end
    close(io)
end

function getgdfprefix(dfnfile::String)
    location = findfirst(".dfn", lowercase(dfnfile))[1] # get file prefix
    dfnfile[1:location-1] # the prefix
end    

function pairinteractionplot(d; varnames=nothing, figsize=(8.5,6), nbins=25, fontsize=8, fbounds=nothing,
        cmap="bone_r", islogpdf=false, showpdf=false, vecofpoints=nothing, vecofpointscolor=nothing,
        scattersize=1, scattercolor="yellowgreen", scatteralpha=1)
    # plot pairs of scatter, d assumed to have realisations in rows
    # vecofpoints should be a vector of vectors if nothing
    @assert !isnothing(varnames)
    nvars = size(d,2)
    f = figure(figsize=figsize)
    T = islogpdf ? x->log(x) : x->x
    for i=1:size(d,2)
        for j=1:i
            c = getrowwise(i,j,nvars)
            ax = subplot(nvars, nvars, c)
            if i==j 
                h = fit(Histogram, d[:,i], nbins=nbins)
                bwidth, bx, denom = getbinsfromhist(h, pdfnormalize=true)
                ax.bar(bx, h.weights./denom, width=bwidth, edgecolor="none", color="yellowgreen")
                ax.set_ylabel("probability")
                ax.yaxis.set_label_position("right")
                ax.tick_params(axis="y", labelright=true, labelleft=false, right=true, left=false)
                if !isnothing(vecofpoints)    
                    for (iv,v) in enumerate(vecofpoints)
                        if isnothing(vecofpointscolor) 
                            ax.plot(v[i]*ones(2),[0, 1], markeredgewidth=3)
                        else
                            ax.plot(v[i]*ones(2),[0, 1], color=vecofpointscolor[iv], linewidth=3)    
                        end
                    end
                end
                if !isnothing(fbounds)    
                    ax.plot(fbounds[i, 1]*ones(2),[0, 1], "--k", linewidth=1)
                    ax.plot(fbounds[i, 2]*ones(2),[0, 1], "--k", linewidth=1)   
                end            
            else
                h = fit(Histogram, (d[:,i], d[:,j]), nbins=nbins)
                if showpdf
                    ax.pcolormesh(h.edges[2], h.edges[1], T.(h.weights/maximum(h.weights)), cmap=cmap)
                else    
                    ax.scatter(d[:,j], d[:,i], s=scattersize, alpha=scatteralpha, c=scattercolor)
                end
                if !isnothing(vecofpoints)
                    for (iv,v) in enumerate(vecofpoints)
                        if isnothing(vecofpointscolor) 
                            ax.plot(v[j], v[i], "+", markersize=10*scattersize, markeredgewidth=3)
                        else
                            ax.plot(v[j], v[i], "+", color=vecofpointscolor[iv], 
                            markersize=10*scattersize, markeredgewidth=3)
                        end        
                    end
                end
                if !isnothing(fbounds)    
                    ax.plot(fbounds[j,1]*ones(2),fbounds[i,:], "--k", linewidth=1)
                    ax.plot(fbounds[j,2]*ones(2),fbounds[i,:], "--k", linewidth=1)
                    ax.plot(fbounds[j,:],fbounds[i,1]*ones(2), "--k", linewidth=1)
                    ax.plot(fbounds[j,:],fbounds[i,2]*ones(2), "--k", linewidth=1)
                end         
                j == 1 && ax.set_ylabel(varnames[i])
                j == 1 || ax.tick_params(labelleft=false)
            end
            (i == nvars) && ax.set_xlabel(varnames[j])
            ax.ticklabel_format(useOffset=false) 
            (i == nvars) ||  ax.tick_params(labelbottom=false)
        end
    end
    ax = f.axes
    for i = size(d,2)-1:-1:1
    # align x axes    
        for j = i:-1:1
            this = firstval(i) + j-1
            next = firstval(i+1) + j-1
            ax[this].sharex(ax[next])
        end
    end
    for i = size(d,2):-1:1
        # align y axes, abundance of caution for
        # weird histograms probably won't need this    
        for j = i-2:-1:1
            this = firstval(i) + j-1
            next = this + 1
            ax[this].sharey(ax[next])
        end
    end    
    # very last histogram shares x 
    # with y axis of previous
    ax[end].set_xlim(ax[end-1].get_ylim())
    nicenup(f, fsize=fontsize)
end

getrowwise(i,j,nvars) = (i-1)*nvars+j
getcolwise(i,j,nvars) = (j-1)*nvars+i
firstval(n) = n == 1 ? 1 : firstval(n-1) + n-1 # index number when filling rowwise upto diagonals
end # module CommonToAll