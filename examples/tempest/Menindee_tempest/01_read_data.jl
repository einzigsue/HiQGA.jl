using PyPlot, Revise, transD_GP
soundings = transD_GP.TEMPEST1DInversion.read_survey_files(
	fname_dat= "TEMPEST_25Hz_Laura_Menindee_EM_Survey_altitude_lines_lidar_elevation.dat",
	fname_specs_halt="electronics_halt.jl",
	frame_height = 19,
	frame_dz = 31,
	frame_dx = 29,
	frame_dy = 30,
	Hxs = [36, 50],
	Hzs = [75, 89],
	Hxp = 69,
	Hzp = 108,
	units = 1e-15,
	yaw_rx = 34,
	pitch_rx = 32,
	roll_rx = 33,
	yaw_tx = 22,
	pitch_tx = 20,
	roll_tx = 21,
	X = 12,
	Y = 13,
	Z = 16,
	fid = 3,
	linenum = 1,
	makesounding=true)
