using Revise
using PyPlot, DelimitedFiles, Random, Statistics, transD_GP

Random.seed!(23)

zfixed = [-1e5]
ρfixed = [1e12]
nmax = 100

zstart = 0.0
extendfrac, dz = 1.06, 2.
zall, znall, zboundaries = transD_GP.setupz(zstart, extendfrac, dz=dz, n=40)
z, ρ, nfixed = transD_GP.makezρ(zboundaries; zfixed=zfixed, ρfixed=ρfixed)
@info z

## geometry parameters for tempest
zTx = -120
zRx = -82
x_rx = -113.0
y_rx = 0.
rx_roll = 0.
rx_pitch = 2.
rx_yaw = 0.
tx_roll = 0.
tx_pitch = 0.
tx_yaw = 0.
# electronics and stuff
include("electronics_halt.jl")
## fill in detail in ohm-m
ρ[(z.>=zstart) .& (z.<50)] .= 20.
ρ[(z.>=50) .&(z.<80)] .= 10.
ρ[(z.>=80) .&(z.<100)] .= 20
ρ[(z.>=100) .&(z.<200)] .= 50
ρ[(z.>=200) .&(z.<250)] .= 80
ρ[(z.>=250)]            .= 150

## create total field operator (required for nuisance inversion)
tempest = transD_GP.TEMPEST1DInversion.Bfield(
    zTx = zTx, zRx = zRx, x_rx = x_rx, y_rx = y_rx,
    rx_roll = rx_roll, rx_pitch = rx_pitch, rx_yaw = rx_yaw,
    tx_roll = tx_roll, tx_pitch = tx_pitch, tx_yaw = tx_yaw,
	ramp = ramp, times = times,
	z=z,
	ρ=ρ,
	addprimary = false)

# plot before adding noise
transD_GP.TEMPEST1DInversion.plotmodelfield!(tempest,z,ρ)
## compute noisy data to invert
# remember noise in electronics_halt.jl are in fT!!
transD_GP.TEMPEST1DInversion.set_noisy_data!(tempest, z, ρ,
                        noisefracx = 0.02, noisefracz = 0.02,
                        halt_X = Hx_add_noise*1e-15, halt_Z = Hz_add_noise*1e-15)
