workspace()

include("file_io/read_key.jl")
include("file_io/write_key.jl")

include("file_io/read_input.jl")
include("file_io/read_output_generation.jl")
include("file_io/write_output.jl")

include("characteristic_matrix_ellipse.jl")
include("characteristic_matrix_ellipsoid.jl")
include("inverse_characteristic_matrix_ellipsoid.jl")
include("rotation_matrix.jl")
include("inverse_rotation_matrix.jl")

include("generate_cell_lists.jl")
include("overlap_cuboid_binary.jl")
include("axis_aligned_bounding_box.jl")
include("intersect_box_box.jl")


foo = @__FILE__
@eval @everywhere f = $foo
#@everywhere println(f)
@everywhere (program_file_dir, program_file_name) = splitdir(f)
@everywhere include(joinpath(program_file_dir, "diffuse.jl"))
@everywhere include(joinpath(program_file_dir, "signed_distance_mod.jl"))
@everywhere include(joinpath(program_file_dir, "position_mod.jl"))


function run_diffusion()
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)

	# Start time.
	t_start_ns::Int64 = convert(Int64, time_ns())

	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)

	# Assert that input is file and store path.
	input_file_path::String = ""
	if isfile(ARGS[1])
		input_file_path = ARGS[1]
	else
		println("No input file specified or specified input file does not exist. Aborting.")
		return nothing
	end

	# Read diffusion input from file.
	println(join(("Reading diffusion input from file ", input_file_path, "...")))

	(	solid_structure_path::String,
		diffusion_map_path::String,
		D0::Float64,
		D1::Float64,
		deltat_coarse::Float64,
		number_of_time_points_coarse::Int64,
		number_of_time_points_fine_per_coarse::Int64,
		number_of_diffusers::Int64,
		number_of_cells_x::Int64,
		number_of_cells_y::Int64,
		number_of_cells_z::Int64,
		boundary_condition::String,
		output_file_path::String) = read_input(input_file_path)

	# Read solid structure from file.
	particle_type_ss::String = ""
	R_ss::Array{Float64, 2} = zeros(0, 0)
	Lx::Float64 = 0.0
	Ly::Float64 = 0.0
	Lz::Float64 = 0.0
	phi::Float64 = 0.0
	X_ss::Array{Float64, 1} = zeros(0)
	Y_ss::Array{Float64, 1} = zeros(0)
	Z_ss::Array{Float64, 1} = zeros(0)
	Q0_ss::Array{Float64, 1} = zeros(0)
	Q1_ss::Array{Float64, 1} = zeros(0)
	Q2_ss::Array{Float64, 1} = zeros(0)
	Q3_ss::Array{Float64, 1} = zeros(0)
	if !isempty(solid_structure_path)
		println(join(("Reading solid structure from file ", solid_structure_path, "...")))
		(	particle_type_ss,
			R_ss,
			Lx,
			Ly,
			Lz,
			phi,
			X_ss,
			Y_ss,
			Z_ss,
			Q0_ss,
			Q1_ss,
			Q2_ss,
			Q3_ss,
			execution_time_generation_ss::Float64) = read_output_generation(solid_structure_path)
	end
	number_of_particles_ss::Int64 = length(X_ss)

	# Read diffusion map from file.
	particle_type_dm::String = ""
	R_dm::Array{Float64, 2} = zeros(0, 0)
	X_dm::Array{Float64, 1} = zeros(0)
	Y_dm::Array{Float64, 1} = zeros(0)
	Z_dm::Array{Float64, 1} = zeros(0)
	Q0_dm::Array{Float64, 1} = zeros(0)
	Q1_dm::Array{Float64, 1} = zeros(0)
	Q2_dm::Array{Float64, 1} = zeros(0)
	Q3_dm::Array{Float64, 1} = zeros(0)
	if !isempty(diffusion_map_path)
		println(join(("Reading solid structure from file ", diffusion_map_path, "...")))
		(	particle_type_dm,
			R_dm,
			Lx,
			Ly,
			Lz,
			phi,
			X_dm,
			Y_dm,
			Z_dm,
			Q0_dm,
			Q1_dm,
			Q2_dm,
			Q3_dm,
			execution_time_generation_dm::Float64) = read_output_generation(diffusion_map_path)
	end
	number_of_particles_dm::Int64 = length(X_dm)

	# Characteristic/rotation matrix entries.
	a11::Float64 = 0.0
	a12::Float64 = 0.0
	a13::Float64 = 0.0
	a21::Float64 = 0.0
	a22::Float64 = 0.0
	a23::Float64 = 0.0
	a31::Float64 = 0.0
	a32::Float64 = 0.0
	a33::Float64 = 0.0

	A11_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	A12_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	A13_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	A21_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	A22_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	A23_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	A31_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	A32_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	A33_ss::Array{Float64, 1} = zeros(number_of_particles_ss)
	if particle_type_ss == "ellipse"
		for current_particle = 1:number_of_particles_ss
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_rotation_matrix(Q0_ss[current_particle], Q1_ss[current_particle], Q2_ss[current_particle], Q3_ss[current_particle])
			A11_ss[current_particle] = a11
			A12_ss[current_particle] = a12
			A13_ss[current_particle] = a13
			A21_ss[current_particle] = a21
			A22_ss[current_particle] = a22
			A23_ss[current_particle] = a23
			A31_ss[current_particle] = a31
			A32_ss[current_particle] = a32
			A33_ss[current_particle] = a33
		end
	elseif particle_type_ss == "ellipsoid"
		for current_particle = 1:number_of_particles_ss
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_characteristic_matrix_ellipsoid(Q0_ss[current_particle], Q1_ss[current_particle], Q2_ss[current_particle], Q3_ss[current_particle], R_ss[current_particle, 1], R_ss[current_particle, 2], R_ss[current_particle, 3])
			A11_ss[current_particle] = a11
			A12_ss[current_particle] = a12
			A13_ss[current_particle] = a13
			A21_ss[current_particle] = a21
			A22_ss[current_particle] = a22
			A23_ss[current_particle] = a23
			A31_ss[current_particle] = a31
			A32_ss[current_particle] = a32
			A33_ss[current_particle] = a33
		end
	elseif particle_type_ss == "cuboid"
		for current_particle = 1:number_of_particles_ss
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_rotation_matrix(Q0_ss[current_particle], Q1_ss[current_particle], Q2_ss[current_particle], Q3_ss[current_particle])
			A11_ss[current_particle] = a11
			A12_ss[current_particle] = a12
			A13_ss[current_particle] = a13
			A21_ss[current_particle] = a21
			A22_ss[current_particle] = a22
			A23_ss[current_particle] = a23
			A31_ss[current_particle] = a31
			A32_ss[current_particle] = a32
			A33_ss[current_particle] = a33
		end
	elseif particle_type_ss == "superellipsoid"
		for current_particle = 1:number_of_particles_ss
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_rotation_matrix(Q0_ss[current_particle], Q1_ss[current_particle], Q2_ss[current_particle], Q3_ss[current_particle])
			A11_ss[current_particle] = a11
			A12_ss[current_particle] = a12
			A13_ss[current_particle] = a13
			A21_ss[current_particle] = a21
			A22_ss[current_particle] = a22
			A23_ss[current_particle] = a23
			A31_ss[current_particle] = a31
			A32_ss[current_particle] = a32
			A33_ss[current_particle] = a33
		end
	end

	A11_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	A12_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	A13_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	A21_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	A22_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	A23_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	A31_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	A32_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	A33_dm::Array{Float64, 1} = zeros(number_of_particles_dm)
	if particle_type_dm == "ellipse"
		for current_particle = 1:number_of_particles_dm
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_rotation_matrix(Q0_dm[current_particle], Q1_dm[current_particle], Q2_dm[current_particle], Q3_dm[current_particle])
			A11_dm[current_particle] = a11
			A12_dm[current_particle] = a12
			A13_dm[current_particle] = a13
			A21_dm[current_particle] = a21
			A22_dm[current_particle] = a22
			A23_dm[current_particle] = a23
			A31_dm[current_particle] = a31
			A32_dm[current_particle] = a32
			A33_dm[current_particle] = a33
		end
	elseif particle_type_dm == "ellipsoid"
		for current_particle = 1:number_of_particles_dm
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_characteristic_matrix_ellipsoid(Q0_dm[current_particle], Q1_dm[current_particle], Q2_dm[current_particle], Q3_dm[current_particle], R_dm[current_particle, 1], R_dm[current_particle, 2], R_dm[current_particle, 3])
			A11_dm[current_particle] = a11
			A12_dm[current_particle] = a12
			A13_dm[current_particle] = a13
			A21_dm[current_particle] = a21
			A22_dm[current_particle] = a22
			A23_dm[current_particle] = a23
			A31_dm[current_particle] = a31
			A32_dm[current_particle] = a32
			A33_dm[current_particle] = a33
		end
	elseif particle_type_dm == "cuboid"
		for current_particle = 1:number_of_particles_dm
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_rotation_matrix(Q0_dm[current_particle], Q1_dm[current_particle], Q2_dm[current_particle], Q3_dm[current_particle])
			A11_dm[current_particle] = a11
			A12_dm[current_particle] = a12
			A13_dm[current_particle] = a13
			A21_dm[current_particle] = a21
			A22_dm[current_particle] = a22
			A23_dm[current_particle] = a23
			A31_dm[current_particle] = a31
			A32_dm[current_particle] = a32
			A33_dm[current_particle] = a33
		end
	elseif particle_type_dm == "superellipsoid"
		for current_particle = 1:number_of_particles_dm
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_rotation_matrix(Q0_dm[current_particle], Q1_dm[current_particle], Q2_dm[current_particle], Q3_dm[current_particle])
			A11_dm[current_particle] = a11
			A12_dm[current_particle] = a12
			A13_dm[current_particle] = a13
			A21_dm[current_particle] = a21
			A22_dm[current_particle] = a22
			A23_dm[current_particle] = a23
			A31_dm[current_particle] = a31
			A32_dm[current_particle] = a32
			A33_dm[current_particle] = a33
		end
	end

	# Create cell lists.
	deltat_fine::Float64 = deltat_coarse / convert(Float64, number_of_time_points_fine_per_coarse)
	sigma::Float64 = sqrt(2.0 * max(D0, D1) * deltat_fine)
	cell_overlap::Float64 = 6.0 * sigma
	cell_lists_ss::Array{Array{Int64, 1}, 3} = generate_cell_lists(	particle_type_ss,
																	R_ss,
																	Lx,
																	Ly,
																	Lz,
																	X_ss,
																	Y_ss,
																	Z_ss,
																	Q0_ss,
																	Q1_ss,
																	Q2_ss,
																	Q3_ss,
																	number_of_cells_x,
																	number_of_cells_y,
																	number_of_cells_z,
																	cell_overlap)
	cell_lists_dm::Array{Array{Int64, 1}, 3} = generate_cell_lists(	particle_type_dm,
																	R_dm,
																	Lx,
																	Ly,
																	Lz,
																	X_dm,
																	Y_dm,
																	Z_dm,
																	Q0_dm,
																	Q1_dm,
																	Q2_dm,
																	Q3_dm,
																	number_of_cells_x,
																	number_of_cells_y,
																	number_of_cells_z,
																	cell_overlap)

	# Simulate diffusion.
	number_of_workers::Int64 = nworkers() # This is determined by the the '-p' input flag to Julia.
	number_of_diffusers_per_worker::Array{Int64, 1} = convert(Array{Int64, 1}, floor(number_of_diffusers / number_of_workers) * ones(number_of_workers))
	number_of_diffusers_remaining::Int64 = number_of_diffusers - sum(number_of_diffusers_per_worker)
	number_of_diffusers_per_worker[1:number_of_diffusers_remaining] += 1
	println(number_of_diffusers_per_worker)

	output::Array{Float64, 1} = zeros(4 * number_of_time_points_coarse + 1)
	output = @parallel (+) for current_worker = 1:number_of_workers
		diffuse(	particle_type_ss, R_ss, X_ss, Y_ss, Z_ss, Q0_ss, Q1_ss, Q2_ss, Q3_ss, A11_ss, A12_ss, A13_ss, A21_ss, A22_ss, A23_ss, A31_ss, A32_ss, A33_ss, cell_lists_ss,
					particle_type_dm, R_dm, X_dm, Y_dm, Z_dm, Q0_dm, Q1_dm, Q2_dm, Q3_dm, A11_dm, A12_dm, A13_dm, A21_dm, A22_dm, A23_dm, A31_dm, A32_dm, A33_dm, cell_lists_dm,
					Lx,	Ly,	Lz,	D0, D1, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse,
					number_of_diffusers_per_worker[current_worker], boundary_condition)
	end

	# Process output.
	msd::Array{Float64, 1} = output[1:number_of_time_points_coarse] ./ convert(Float64, 3 * number_of_diffusers)
	msd_x::Array{Float64, 1} = output[number_of_time_points_coarse+1:2*number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	msd_y::Array{Float64, 1} = output[2*number_of_time_points_coarse+1:3*number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	msd_z::Array{Float64, 1} = output[3*number_of_time_points_coarse+1:4*number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	D_empirical::Float64 = output[end] / (3.0 * convert(Float64, number_of_diffusers * (number_of_time_points_coarse-1) * number_of_time_points_fine_per_coarse) * 2.0 * deltat_fine)
	t_finish_ns::Int64 = convert(Int64, time_ns())
	t_exec::Float64 = convert(Float64, t_finish_ns - t_start_ns) / 1e9

	# Write output.
	write_output(
		output_file_path,
		D_empirical,
		deltat_coarse,
		number_of_time_points_coarse,
		msd,
		msd_x,
		msd_y,
		msd_z,
		t_exec)
	println(join(("Output written to ", output_file_path, ".")))
	println("Finished.")
	#println(msd_z[end]/(2.0*(convert(Float64, number_of_time_points_coarse-1) * deltat_coarse)))

	nothing
end

run_diffusion()
