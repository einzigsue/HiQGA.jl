# wrapper to integrate SMRPI 1D forward model with
# HiQGA/transD_GP
module SMRPI
using ..AbstractOperator, ..CommonToAll
import ..AbstractOperator.get_misfit
import ..Model, ..Options
using SNMRForward, Random

mutable struct SMRSounding <: Operator1D
    V0 :: Vector{<:Real} #sounding curve
    ϕ :: Vector{<:Real} #phases
    σ_V0 :: Vector{<:Real}
    σ_ϕ :: Vector{<:Real}
    Fm :: SNMRForward.MRSForward
end

function newSMRSounding(V0, ϕ, σ_V0, σ_ϕ, Fm)
    (length(V0) != length(ϕ) ||
    length(V0) != length(σ_V0) ||
    length(ϕ) != length(σ_ϕ)) && 
    throw(ArgumentError("V0, ϕ and associated errors must have same length"))

    SMRSounding(V0, ϕ, σ_V0, σ_ϕ, Fm)
end

function get_misfit(m::Model, opt::Options, S::SMRSounding)
    opt.debug && return 0.0
    get_misfit(10 .^m.fstar[:], S)
end
# above defined function and type signature MUST be defined

function get_misfit(w::Vector{<:Real}, S::SMRSounding)
    response = SNMRForward.forward(S.Fm, w)
    Vres = abs.(response)
    ϕres = angle.(response)
    residual = [(S.V0 .- Vres)./S.σ_V0; (S.ϕ .- ϕres)./S.σ_ϕ]
    residual' * residual / 2
end

function create_synthetic(w::Vector{<:Real}, σ::Vector{<:Real}, t::Vector{<:Real},
            Be::Real, ϕ::Real, R::Real, zgrid::Vector{<:Real}, qgrid::Vector{<:Real}
    ; noise_frac = 0.05)
    ct = SNMRForward.ConductivityModel(σ, t)

    F = SNMRForward.MRSForward(R, zgrid, qgrid, ϕ, Be, ct)

    synth_data = SNMRForward.forward(F,w)
    σ_V0 = noise_frac * abs.(synth_data)
    σ_ϕ = noise_frac * ones(size(synth_data))
    noisy_V0 = abs.(synth_data) .+ σ_V0 .* randn(size(abs.(synth_data)))
    noisy_ϕ = angle.(synth_data) .+ σ_ϕ .* randn(size(abs.(synth_data)))

    newSMRSounding(noisy_V0, noisy_ϕ, σ_V0, σ_ϕ, F)
end


end