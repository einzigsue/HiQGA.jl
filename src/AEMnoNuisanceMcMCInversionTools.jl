module AEMnoNuisanceMcMCInversionTools
using ..AbstractOperator, ..CommonToAll
import ..AbstractOperator.makeoperator
import ..AbstractOperator.loopacrossAEMsoundings
import ..AbstractOperator.plotmodelfield!
import ..AbstractOperator.getndata
import ..Options, ..OptionsStat
export makeAEMoperatorandoptions, loopacrossAEMsoundings, summaryAEMimages, plotindividualAEMsoundings
import ..main # McMC function
using ..SoundingDistributor
using Distributed, Dates, Statistics, DelimitedFiles, PyPlot, Random
# plotting stuff
function summaryAEMimages(soundings::Array{S, 1}, opt::Options;
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        dz = -1,
                        dr = 10,
                        zall=[-1.],
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap="turbo",
                        figsize=(6,10),
                        topowidth=2,
                        idx = nothing,
                        omitconvergence = false,
                        useML = false,
                        preferEright = false,
                        preferNright = false,
                        saveplot = true,
                        yl = nothing,
                        showplot = true,
                        showmean = false,
                        dpi = 300) where S<:Sounding
    
    linestartidx = splitsoundingsbyline(soundings)                    
    nlines = length(linestartidx)                   
    for i in 1:nlines
        a = linestartidx[i]
        b = i != nlines ?  linestartidx[i+1]-1 : length(soundings)
        summaryimages(soundings[a:b], opt; qp1, qp2, burninfrac, zall,dz, dr, 
            fontsize, vmin, vmax, cmap, figsize, topowidth, idx=idx, omitconvergence, useML, 
            preferEright, showplot, preferNright, saveplot, yl, dpi, showmean)
    end
    nothing    
end

function summaryimages(soundings::Array{S, 1}, opt::Options;
                        qp1=0.05,
                        qp2=0.95,
                        burninfrac=0.5,
                        dz = -1,
                        dr = 10,
                        zall = [-1.],
                        fontsize = 10,
                        vmin = -2,
                        vmax = 0.5,
                        cmap="turbo",
                        figsize=(6,10),
                        topowidth=2,
                        idx = nothing,
                        omitconvergence = false,
                        useML = false,
                        preferEright = false,
                        preferNright = false,
                        saveplot = false,
                        yl = nothing,
                        showplot = true,
                        showmean = false,
                        dpi = 300) where S<:Sounding
    @assert !(preferNright && preferEright) # can't prefer both labels to the right
    pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd  = summarypost(soundings, opt; zall, qp1, qp2, burninfrac, useML)

    phgrid, plgrid, pmgrid, σmeangrid, ∇zmeangrid,
    ∇zsdgrid, gridx, gridz, topofine, R, Z = makesummarygrid(soundings, pl, pm, ph, ρmean,
                                                            vdmean, vddev, zall, dz, dr=dr)

    lname = "Line $(soundings[1].linenum)"
    plotsummarygrids1(soundings, σmeangrid, phgrid, plgrid, pmgrid, gridx, gridz, topofine, R, Z, χ²mean, χ²sd, lname; qp1, qp2,
                        figsize, fontsize, cmap, vmin, vmax, 
                        topowidth, idx, omitconvergence, useML,
                        preferEright, preferNright, saveplot, showplot, dpi,
                        yl, showmean)                  
end

function summarypost(soundings::Vector{S}, opt::Options;
            qp1=0.05,
            qp2=0.95,
            burninfrac=0.5,
            zall = [-1.],
            useML=false) where S<:Sounding

    @assert length(zall) != 1

    linename = "_line_$(soundings[1].linenum)_summary.txt"
    fnames = ["rho_low", "rho_mid", "rho_hi", "rho_avg",
              "ddz_mean", "ddz_sdev", "phid_mean", "phid_sdev"].*linename
    if isfile(fnames[1])
        @warn fnames[1]*" exists, reading stored values"
        pl, pm, ph, ρmean,
        vdmean, vddev, χ²mean, χ²sd, = map(x->readdlm(x), fnames)
    else
        # this is a dummy operator for plotting
        aem, = makeoperator(soundings[1])
        # now get the posterior marginals
        opt.xall[:] .= zall
        pl, pm, ph, ρmean, vdmean, vddev = map(x->zeros(length(zall), length(soundings)), 1:6)
        χ²mean, χ²sd = zeros(length(soundings)), zeros(length(soundings))
        for idx = 1:length(soundings)
            @info "$idx out of $(length(soundings))\n"
            opt.fdataname = soundings[idx].sounding_string*"_"
            opt.xall[:] .= zall
            pl[:,idx], pm[:,idx], ph[:,idx], ρmean[:,idx],
            vdmean[:,idx], vddev[:,idx] = CommonToAll.plot_posterior(aem, opt, burninfrac=burninfrac,
                                                    qp1=qp1, qp2=qp2,
                                                    doplot=false)
            χ² = 2*CommonToAll.assembleTat1(opt, :U, temperaturenum=1, burninfrac=burninfrac)
            ndata = getndata(soundings[idx])
            χ²mean[idx] = mean(χ²)/ndata
            χ²sd[idx]   = std(χ²)/ndata
            if useML
                # this is approximate as HM and LM have different ML factors sampled 
                χ²mean[idx] = exp(χ²mean[idx]-log(ndata))
                χ²sd[idx]   = exp(χ²sd[idx]-log(ndata)) # I think, need to check
            end    
        end
        # write in grid format
        for (fname, vals) in Dict(zip(fnames, [pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd]))
            writedlm(fname, vals)
        end
        # write in x, y, z, rho format
        for (i, d) in enumerate([pl, pm, ph, ρmean])
            xyzrho = makearray(soundings, d, zall)
            writedlm(fnames[i][1:end-4]*"_xyzrho.txt", xyzrho)
        end
    end
    pl, pm, ph, ρmean, vdmean, vddev, χ²mean, χ²sd
end

function getndata(d)
    select = .!isnan.(d)
    ndata  = sum(select)
    ndata, select
end

function plotindividualAEMsoundings(soundings::Vector{S}, aem_in::Operator1D, opt_in::Options, idxplot::Vector{Int};
    zall=[-1.],
    burninfrac=0.5,
    nbins = 50,
    figsize  = (6,6),
    omittemp = true,
    showslope = false,
    plotmean = false,
    pdfclim = nothing,
    qp1=0.05,
    qp2=0.95,
    rseed = 123,
    computeforwards = false,
    nforwards = 100) where S<:Sounding
    
    @assert length(zall) != 1
    opt = deepcopy(opt_in)
    opt.xall[:] = zall
    for idx = 1:length(soundings)
        if in(idx, idxplot)
            @info "Sounding number: $idx"
            aem = makeoperator(aem_in, soundings[idx])
            opt.fdataname = soundings[idx].sounding_string*"_"
            getchi2forall(opt, alpha=0.8, omittemp=omittemp)
            CommonToAll.getstats(opt)
            plot_posterior(aem, opt; burninfrac, nbins, figsize, 
                    showslope, pdfclim, plotmean, qp1, qp2)
            ax = gcf().axes
            ax[1].invert_xaxis()
            if computeforwards
                M = assembleTat1(opt, :fstar, temperaturenum=1, burninfrac=burninfrac)
                Random.seed!(rseed)
                plotmodelfield!(aem, M[randperm(length(M))[1:nforwards]])
            end            
        end
    end
end

# stuff needed for McMC driver code
function make_tdgp_opt(;
                    rseed = nothing,
                    znall = znall,
                    fileprefix = "sounding",
                    nmin = 2,
                    nmax = 40,
                    K = GP.OrstUhn(),
                    demean = true,
                    sdpos = 0.05,
                    sdprop = 0.05,
                    sddc = 0.008,
                    sampledc = false,
                    fbounds = [-0.5 2.5],
                    λ = [2],
                    δ = 0.1,
                    pnorm = 2,
                    save_freq = 25,
                    restart = false
                    )
    sdev_pos = [sdpos*abs(diff([extrema(znall)...])[1])]
    sdev_prop = sdprop*diff(fbounds, dims=2)[:]
    sdev_dc = sddc*diff(fbounds, dims=2)[:]
    xall = permutedims(collect(znall))
    xbounds = permutedims([extrema(znall)...])

    history_mode = "w"
	restart && (history_mode = "a")

    if rseed != nothing
        Random.seed!(rseed)
    end

    opt = OptionsStat(;fdataname = fileprefix*"_",
        nmin, nmax, xbounds, fbounds, xall, λ, δ,
        demean, sdev_prop, sdev_pos, sdev_dc,
        sampledc, pnorm, quasimultid = false, K,
        save_freq, needλ²fromlog=false, updatenonstat=false,
        dispstatstoscreen = false, history_mode)
    opt
end

function makeAEMoperatorandoptions(sounding::Sounding;
                        nmin = 2,
                        nmax = 40,
                        K = GP.OrstUhn(),
                        demean = false,
                        sdpos = 0.05,
                        sdprop = 0.05,
                        sddc = 0.008,
                        sampledc = false,
                        fbounds = [-0.5 2.5],
                        λ = [2],
                        δ = 0.1,
                        save_freq = 25,
                        restart = false,
                        zfixed   = [-1e5],
                        ρfixed   = [1e12],
                        zstart = 0.0,
                        extendfrac = 1.06,
                        dz = 2.,
                        ρbg = 10,
                        nlayers = 40,
                        ntimesperdecade = 10,
                        nfreqsperdecade = 5,
                        showgeomplot = false,
                        plotfield = false,
                        useML = false,
                        modelprimary = false
                        )
    aem, zall, znall,  = makeoperator(sounding;
                        zfixed, ρfixed, zstart, extendfrac,
                        dz = dz, ρbg, nlayers, ntimesperdecade,
                        nfreqsperdecade, useML, showgeomplot,
                        modelprimary, plotfield)

    opt = make_tdgp_opt(znall = znall,
                        fileprefix = sounding.sounding_string,
                        nmin = nmin,
                        nmax = nmax,
                        K = K,
                        demean = demean,
                        sdpos = sdpos,
                        sdprop = sdprop,
                        sddc = sddc,
                        sampledc = sampledc,
                        fbounds = fbounds,
                        save_freq = save_freq,
                        λ = λ,
                        δ = δ,
                        restart = restart
                        )
    aem, opt, zall
end

# Driver code for McMC inversion with no nuisances
function loopacrossAEMsoundings(soundings::Array{S, 1}, aem_in::Operator1D, opt_in::Options;
                            Tmax               = -1,
                            nsamples           = -1,
                            nchainsatone       =  1,
                            nchainspersounding = -1,
                            ppn                = -1) where S<:Sounding

    @assert ppn != -1
    @assert nchainspersounding != -1
    @assert nsamples != - 1
    @assert Tmax != -1

    nsoundings = length(soundings)
    nsequentialiters, nparallelsoundings = splittasks(soundings; nchainspersounding, ppn)
    opt = deepcopy(opt_in)
    
    for iter = 1:nsequentialiters
        ss = getss(iter, nsequentialiters, nparallelsoundings, nsoundings)
        @info "soundings in loop $iter of $nsequentialiters", ss
        @sync for (i, s) in Iterators.reverse(enumerate(ss))
            pids = getpids(i, nchainspersounding)
            @info "pids in sounding $s:", pids
            
            aem = makeoperator(aem_in, soundings[s])
            opt.fdataname = soundings[s].sounding_string*"_"

            @async remotecall_wait(main, pids[1], opt, aem, collect(pids[2:end]),
                                    Tmax         = Tmax,
                                    nsamples     = nsamples,
                                    nchainsatone = nchainsatone)

        end # @sync
        @info "done $iter out of $nsequentialiters at $(Dates.now())"
    end
end

end