## read the soundings
using PyPlot, Revise, transD_GP
soundings = transD_GP.TEMPEST1DInversion.read_survey_files(
	fname_dat= "91010901.asc",
	fname_specs_halt="electronics_halt.jl",
	frame_height = 20,
	frame_dz = 32,
	frame_dx = 30,
	frame_dy = 31,
	Hxs = [37, 51],
	Hzs = [76, 90],
	Hxp = 70,
	Hzp = 109,
	units = 1e-15,
	yaw_rx = 35,
	pitch_rx = 33,
	roll_rx = 34,
	yaw_tx = 23,
	pitch_tx = 21,
	roll_tx = 22,
	X = 12,
	Y = 13,
	Z = 18,
	fid = 3,
	startfrom = 1,
	skipevery = 2,
	linenum = 1,
	makesounding=true)

