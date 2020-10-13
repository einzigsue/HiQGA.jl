module SkyTEM1DInversion
import AbstractOperator.get_misfit
using AbstractOperator, AEM_VMD_HMD
using TransD_GP, PyPlot, LinearAlgebra, CommonToAll, MAT, Random

export dBzdt, plotmodelfield!, addnoise_skytem, plotmodelfield!

const μ₀ = 4*pi*1e-7

mutable struct dBzdt<:Operator1D
    dlow       :: Array{Float64, 1}
    dhigh      :: Array{Float64, 1}
    useML      :: Bool
    σlow       :: Array{Float64, 1}
    σhigh      :: Array{Float64, 1}
    Flow       :: AEM_VMD_HMD.HField
    Fhigh      :: AEM_VMD_HMD.HField
    z          :: Array{Float64, 1}
    nfixed     :: Int
    ρ          :: Array{Float64, 1}
    selectlow  :: Array{Bool, 1}
    selecthigh :: Array{Bool, 1}
    ndatalow   :: Int
    ndatahigh  :: Int
end

function dBzdt(Flow       :: AEM_VMD_HMD.HField,
                      Fhigh      :: AEM_VMD_HMD.HField,
                      dlow       :: Array{Float64, 1},
                      dhigh      :: Array{Float64, 1},
                      σlow       :: Array{Float64, 1},
                      σhigh      :: Array{Float64, 1};
                      useML  = false,
                      z = [-1.],
                      ρ = [-1.],
                      nfixed = 1
                    )
    @assert length(Flow.thickness) >= length(z)
    @assert length(Fhigh.thickness) >= length(z)
    @assert size(σlow)  == size(dlow)
    @assert size(σhigh) == size(dhigh)
    ndatalow   = sum(.!isnan.(dlow))
    ndatahigh  = sum(.!isnan.(dhigh))
    selectlow  = .!isnan.(dlow)
    selecthigh = .!isnan.(dhigh)
    dBzdt(dlow, dhigh, useML, σlow, σhigh,
    Flow, Fhigh, z, nfixed, copy(ρ), selectlow, selecthigh, ndatalow, ndatahigh)
end

function getfield!(m::TransD_GP.Model, aem::dBzdt)
    copyto!(aem.ρ, aem.nfixed+1:length(aem.ρ), 10 .^m.fstar, 1:length(m.fstar))
    AEM_VMD_HMD.getfieldTD!(aem.Flow,  aem.z, aem.ρ)
    AEM_VMD_HMD.getfieldTD!(aem.Fhigh, aem.z, aem.ρ)
    nothing
end

function getfield!(m::Array{Float64}, aem::dBzdt)
    copyto!(aem.ρ, aem.nfixed+1:length(aem.ρ), 10 .^m, 1:length(m))
    AEM_VMD_HMD.getfieldTD!(aem.Flow,  aem.z, aem.ρ)
    AEM_VMD_HMD.getfieldTD!(aem.Fhigh, aem.z, aem.ρ)
    nothing
end

function get_misfit(m::TransD_GP.Model, opt::TransD_GP.Options, aem::dBzdt)
    chi2by2 = 0.0
    if !opt.debug
        getfield!(m, aem)
        chi2by2 += getchi2by2(aem.Flow.dBzdt, aem.dlow,
                    aem.σlow, aem.selectlow, aem.useML, aem.ndatalow)
        chi2by2 += getchi2by2(aem.Fhigh.dBzdt, aem.dhigh,
                    aem.σhigh, aem.selecthigh, aem.useML, aem.ndatahigh)
    end
    return chi2by2
end

function getchi2by2(dBzdt, d, σ, select, useML, ndata)
    r, d, s, idx = dBzdt, d, σ, select
    r .= (r - d)./s
    if useML
        chi2by2 = 0.5*ndata[ifreq]*log(norm(r[idx])^2)
    else
        chi2by2 = 0.5*norm(r[idx])^2
    end
end

function plotmodelfield_skytem!(Flow::AEM_VMD_HMD.HField, Fhigh::AEM_VMD_HMD.HField,
                         z::Array{Float64, 1}, ρ::Array{Float64, 1}
                        ;figsize=(10,5))
    f, ax = plt.subplots(1, 2, figsize=figsize)
    ax[1].step(log10.(ρ[2:end]), z[2:end])
    AEM_VMD_HMD.getfieldTD!(Flow, z, ρ)
    AEM_VMD_HMD.getfieldTD!(Fhigh, z, ρ)
    ax[2].loglog(Flow.times,μ₀*Flow.dBzdt, label="low moment")
    ax[2].loglog(Fhigh.times,μ₀*Fhigh.dBzdt, label="high moment")
    ax[1].grid()
    ax[1].invert_yaxis()
    ax[2].grid()
    nicenup(f)
end

function addnoise_skytem(Flow::AEM_VMD_HMD.HField, Fhigh::AEM_VMD_HMD.HField,
                  z::Array{Float64, 1}, ρ::Array{Float64, 1};
                  noisefrac  = 0.05,
                  noisefloorlow = μ₀*1e-14,
                  noisefloorhigh = μ₀*1e-14,
                  dz = -1.,
                  extendfrac = -1.,
                  nfixed = -1,
                  figsize=(12,4)
                  )
    @assert all((nfixed, dz, extendfrac) .> 0)
    AEM_VMD_HMD.getfieldTD!(Flow, z, ρ)
    AEM_VMD_HMD.getfieldTD!(Fhigh, z, ρ)
    dlow  = Flow.dBzdt + noisefrac*abs.(Flow.dBzdt).*randn(size(Flow.dBzdt))
    dhigh = Fhigh.dBzdt + noisefrac*abs.(Fhigh.dBzdt).*randn(size(Fhigh.dBzdt))
    dlow[abs.(dlow).<noisefloorlow] .= NaN
    dhigh[abs.(dhigh).<noisefloorhigh] .= NaN
    σlow = noisefrac*abs.(dlow)
    σhigh = noisefrac*abs.(dhigh)
    plotmodelfield_skytem!(Flow, Fhigh, z, ρ, dlow, dhigh, σlow, σhigh,
                                figsize=figsize, nfixed=nfixed,
                                dz=dz, extendfrac=extendfrac)
    # returned data is dBzdt not H if there is a μ multiplied
    return μ₀.*(dlow, dhigh, σlow, σhigh)
end

function plotmodelfield!(Flow::AEM_VMD_HMD.HField, Fhigh::AEM_VMD_HMD.HField,
                        z, ρ, dlow, dhigh, σlow, σhigh;
                        figsize=(12,4), nfixed=-1, dz=-1., extendfrac=-1.)
    # expects data and noise in units of H, i.e. B/μ
    @assert all((nfixed, dz, extendfrac) .> 0)
    f, ax = plt.subplots(1, 2, figsize=figsize)
    ax[1].step(log10.(ρ[2:end]), z[2:end])
    if dz > 0
        axn = ax[1].twinx()
        ax[1].get_shared_y_axes().join(ax[1],axn)
        axn.step(log10.(ρ[2:end]), z[2:end])
        yt = ax[1].get_yticks()[ax[1].get_yticks().>=0]
        axn.set_yticks(yt)
        axn.set_ylim(ax[1].get_ylim()[end:-1:1])
        axn.set_yticklabels(string.(Int.(round.(getn.(yt .- z[nfixed+1], dz, extendfrac)))))
    end
    AEM_VMD_HMD.getfieldTD!(Flow, z, ρ)
    AEM_VMD_HMD.getfieldTD!(Fhigh, z, ρ)
    ax[2].loglog(Flow.times,μ₀*Flow.dBzdt, label="low moment")
    ax[2].loglog(Fhigh.times,μ₀*Fhigh.dBzdt, label="high moment")
    ax[2].errorbar(Flow.times, μ₀*vec(dlow), yerr = μ₀*2abs.(vec(σlow)),
                        linestyle="none", marker=".", elinewidth=0, capsize=3)
    ax[2].errorbar(Fhigh.times, μ₀*vec(dhigh), yerr = μ₀*2abs.(vec(σhigh)),
                        linestyle="none", marker=".", elinewidth=0, capsize=3)
    ax[1].grid()
    ax[2].grid()
    nicenup(f)
end

function plotmodelfield!(aem::dBzdt, Ρ::Vector{Array{Float64}};
                        figsize=(8,5), dz=-1., onesigma=true,
                        extendfrac=-1., fsize=10, alpha=0.1)
    @assert all((dz, extendfrac) .> 0)
    sigma = onesigma ? 1.0 : 2.0
    f = figure(figsize=figsize)
    ax = Vector{PyPlot.PyObject}(undef, 3)
    ax[1] = subplot(121)
    ρmin, ρmax = extrema(vcat(Ρ...))
    delρ = ρmax - ρmin
    ax[1].set_xlim(ρmin-0.1delρ,ρmax+0.1delρ)
    nfixed, z = aem.nfixed, aem.z
    ax[1].plot([ρmin-0.1delρ,ρmax+0.1delρ], z[nfixed+1]*[1., 1], color="b")
    ax[2] = subplot(122)
    Flow = aem.Flow
    dlow, σlow = aem.dlow, aem.σlow
    Fhigh = aem.Fhigh
    dhigh, σhigh = aem.dhigh, aem.σhigh
    ax[2].errorbar(Flow.times, μ₀*dlow, yerr = μ₀*sigma*abs.(σlow),
                        linestyle="none", marker=".", elinewidth=0, capsize=3, label="low moment")
    ax[2].errorbar(Fhigh.times, μ₀*dhigh, yerr = μ₀*sigma*abs.(σhigh),
                        linestyle="none", marker=".", elinewidth=0, capsize=3, label="high moment")
    for ρ in Ρ
        getfield!(ρ,  aem)
        Flow.dBzdt[.!aem.selectlow] .= NaN
        Fhigh.dBzdt[.!aem.selecthigh] .= NaN
        ax[1].step(log10.(aem.ρ[2:end]), aem.z[2:end], "-k", alpha=alpha)
        ax[2].loglog(Flow.times,μ₀*Flow.dBzdt, "k", alpha=alpha, markersize=2)
        ax[2].loglog(Fhigh.times,μ₀*Fhigh.dBzdt, "k", alpha=alpha, markersize=2)
    end
    ax[1].grid()
    ax[1].set_ylabel("Depth m")
    ax[1].plot(xlim(), z[nfixed+1]*[1, 1], "--k")
    if dz > 0
        axn = ax[1].twinx()
        ax[1].get_shared_y_axes().join(ax[1],axn)
        yt = ax[1].get_yticks()[ax[1].get_yticks().>=z[nfixed+1]]
        axn.set_yticks(yt)
        axn.set_ylim(ax[1].get_ylim()[end:-1:1])
        axn.set_yticklabels(string.(Int.(round.(getn.(yt .- z[nfixed+1], dz, extendfrac)))))
    end
    axn.set_ylabel("Depth index", rotation=-90)
    ax[1].set_xlabel("Log₁₀ρ")
    ax[1].set_title("Model")
    ax[2].set_ylabel(L"dBzdt \; V/(A.m^4)")
    ax[2].set_xlabel("Time (s)")
    ax[2].set_title("Transient response")
    ax[2].legend()
    ax[2].grid()
    ax[1].invert_xaxis()
    nicenup(f, fsize=fsize)
end

function makeoperator(fdataname::String;
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
                       plotfield = false)
    @assert extendfrac > 1.0
    @assert dz > 0.0
    @assert ρbg > 0.0
    @assert nlayers > 1
    nmax = nlayers+1

    zall, znall, zboundaries = setupz(zstart, extendfrac, dz=dz, n=nlayers, showplot=showgeomplot)
    z, ρ, nfixed = makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
    ρ[z.>=zstart] .= ρbg
    ##  geometry and modeling parameters
    file = matopen(fdataname)
    rRx = read(file, "rRx")
    zRxLM = read(file, "LM_zRx")
    zTxLM = read(file, "LM_zTx")
    zRxHM = read(file, "HM_zRx")
    zTxHM = read(file, "HM_zTx")
    rTx = read(file, "rTxLoop")
    lowpassfcs = read(file, "lowPassFilters")
    # Note that the receiver depth needs to be in same model layer as transmitter.
    ## LM times and ramp
    LM_times = read(file, "LM_times")[:]
    LM_ramp = read(file, "LM_ramp")
    ## HM times and ramp
    HM_times = read(file, "HM_times")[:]
    HM_ramp = read(file, "HM_ramp")
    ## LM operator
    Flm = AEM_VMD_HMD.HFieldDHT(
                          ntimesperdecade = ntimesperdecade,
                          nfreqsperdecade = nfreqsperdecade,
                          lowpassfcs = lowpassfcs,
                          times  = LM_times,
                          ramp   = LM_ramp,
                          nmax   = nmax,
                          zTx    = zTxLM,
                          rRx    = rRx,
                          rTx    = rTx,
                          zRx    = zRxLM)
    ## HM operator
    Fhm = AEM_VMD_HMD.HFieldDHT(
                          ntimesperdecade = ntimesperdecade,
                          nfreqsperdecade = nfreqsperdecade,
                          lowpassfcs = lowpassfcs,
                          times  = HM_times,
                          ramp   = HM_ramp,
                          nmax   = nmax,
                          zTx    = zTxHM,
                          rRx    = rRx,
                          rTx    = rTx,
                          zRx    = zRxHM)
    ## data and high altitude noise
    LM_data = read(file, "d_LM")
    HM_data = read(file, "d_HM")
    LM_noise = read(file, "sd_LM")
    HM_noise = read(file, "sd_HM")
    ## create operator
    dlow, dhigh, σlow, σhigh = (LM_data, HM_data, LM_noise, HM_noise)./μ₀
    aem = dBzdt(Flm, Fhm, vec(dlow), vec(dhigh),
                                      vec(σlow), vec(σhigh), z=z, ρ=ρ, nfixed=nfixed)
    plotfield && plotmodelfield!(Flm, Fhm, z, ρ, dlow, dhigh, σlow, σhigh;
                          figsize=(12,4), nfixed=nfixed, dz=dz, extendfrac=extendfrac)
    aem, znall
end

function make_tdgp_statmode_opt(;
                    rseed = nothing,
                    znall = znall,
                    fileprefix = "sounding",
                    nmin = 2,
                    nmax = 40,
                    K = GP.Mat32(),
                    demean = true,
                    sdpos = 0.05,
                    sdprop = 0.05,
                    fbounds = [-0.5 2.5],
                    λ = [2],
                    δ = 0.1,
                    pnorm = 2,
                    save_freq = 25
                    )
    sdev_pos = [sdpos*abs(diff([extrema(znall)...])[1])]
    sdev_prop = sdprop*diff(fbounds, dims=2)[:]
    xall = permutedims(collect(znall))
    xbounds = permutedims([extrema(znall)...])

    updatenonstat = false
    needλ²fromlog = false
    if rseed != nothing
        Random.seed!(rseed)
    end
    opt = TransD_GP.OptionsStat(fdataname = fileprefix*"_",
                            nmin = nmin,
                            nmax = nmax,
                            xbounds = xbounds,
                            fbounds = fbounds,
                            xall = xall,
                            λ = λ,
                            δ = δ,
                            demean = demean,
                            sdev_prop = sdev_prop,
                            sdev_pos = sdev_pos,
                            pnorm = pnorm,
                            quasimultid = false,
                            K = K,
                            save_freq = save_freq,
                            needλ²fromlog = needλ²fromlog,
                            updatenonstat = updatenonstat
                            )

    ## Initialize options for the dummy nonstationary properties GP
    optdummy = TransD_GP.OptionsNonstat(opt,
                            nmin = 2,
                            nmax = 3,
                            fbounds = fbounds,
                            δ = δ,
                            demean = demean,
                            sdev_prop = sdev_prop,
                            sdev_pos = sdev_pos,
                            pnorm = pnorm,
                            K = K
                            )

    opt, optdummy

end


end
