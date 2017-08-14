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


foo = @__FILE__
@eval @everywhere f = $foo
@everywhere println(f)
@everywhere (program_file_dir, program_file_name) = splitdir(f)
@everywhere include(joinpath(program_file_dir, "diffuse.jl"))
@everywhere include(joinpath(program_file_dir, "signed_distance_mod.jl"))
@everywhere include(joinpath(program_file_dir, "position_mod.jl"))
@everywhere include(joinpath(program_file_dir, "generate_cell_lists.jl"))
@everywhere include(joinpath(program_file_dir, "axis_aligned_bounding_box.jl"))
@everywhere include(joinpath(program_file_dir, "intersect_box_box.jl"))

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

	(	output_generation_path::String,
		D0::Float64,
		deltat_coarse::Float64,
		number_of_time_points_coarse::Int64,
		number_of_time_points_fine_per_coarse::Int64,
		number_of_diffusers::Int64,
		number_of_cells_x::Int64,
		number_of_cells_y::Int64,
		number_of_cells_z::Int64,
		boundary_condition::String,
		output_file_path::String) = read_input(input_file_path)

	# Read generation output from file.
	println(join(("Reading generation output from file ", output_generation_path, "...")))
	(	particle_type::String,
		R::Array{Float64, 2},
		Lx::Float64,
		Ly::Float64,
		Lz::Float64,
		phi::Float64,
		X::Array{Float64, 1},
		Y::Array{Float64, 1},
		Z::Array{Float64, 1},
		Q0::Array{Float64, 1},
		Q1::Array{Float64, 1},
		Q2::Array{Float64, 1},
		Q3::Array{Float64, 1},
		execution_time_diffusion::Float64) = read_output_generation(output_generation_path)
	number_of_particles::Int64 = length(X)

	# Characteristic/rotation matrix entries.
	A11::Array{Float64, 1} = zeros(number_of_particles)
	A12::Array{Float64, 1} = zeros(number_of_particles)
	A13::Array{Float64, 1} = zeros(number_of_particles)
	A21::Array{Float64, 1} = zeros(number_of_particles)
	A22::Array{Float64, 1} = zeros(number_of_particles)
	A23::Array{Float64, 1} = zeros(number_of_particles)
	A31::Array{Float64, 1} = zeros(number_of_particles)
	A32::Array{Float64, 1} = zeros(number_of_particles)
	A33::Array{Float64, 1} = zeros(number_of_particles)

	a11::Float64 = 0.0
	a12::Float64 = 0.0
	a13::Float64 = 0.0
	a21::Float64 = 0.0
	a22::Float64 = 0.0
	a23::Float64 = 0.0
	a31::Float64 = 0.0
	a32::Float64 = 0.0
	a33::Float64 = 0.0
	if particle_type == "ellipse"
		for current_particle = 1:number_of_particles
			#(a11, a12, a13, a21, a22, a23, a31, a32, a33) = characteristic_matrix_ellipse(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle], R[current_particle, 1], R[current_particle, 2])
			#A11[current_particle] = a11
			#A12[current_particle] = a12
			#A13[current_particle] = a13
			#A21[current_particle] = a21
			#A22[current_particle] = a22
			#A23[current_particle] = a23
			#A31[current_particle] = a31
			#A32[current_particle] = a32
			#A33[current_particle] = a33
		end
	elseif particle_type == "ellipsoid"
		for current_particle = 1:number_of_particles
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_characteristic_matrix_ellipsoid(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle], R[current_particle, 1], R[current_particle, 2], R[current_particle, 3])
			A11[current_particle] = a11
			A12[current_particle] = a12
			A13[current_particle] = a13
			A21[current_particle] = a21
			A22[current_particle] = a22
			A23[current_particle] = a23
			A31[current_particle] = a31
			A32[current_particle] = a32
			A33[current_particle] = a33
		end
	elseif particle_type == "cuboid"
		for current_particle = 1:number_of_particles
			(a11, a12, a13, a21, a22, a23, a31, a32, a33) = inverse_rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
			A11[current_particle] = a11
			A12[current_particle] = a12
			A13[current_particle] = a13
			A21[current_particle] = a21
			A22[current_particle] = a22
			A23[current_particle] = a23
			A31[current_particle] = a31
			A32[current_particle] = a32
			A33[current_particle] = a33
		end
	end

	# Create cell lists.
	deltat_fine::Float64 = deltat_coarse / convert(Float64, number_of_time_points_fine_per_coarse)
	sigma::Float64 = sqrt(2.0 * D0 * deltat_fine)
	cell_overlap::Float64 = 6.0 * sigma
	cell_lists::Array{Array{Int64, 1}, 3} = generate_cell_lists(particle_type,
																R,
																Lx,
																Ly,
																Lz,
																X,
																Y,
																Z,
																Q0,
																Q1,
																Q2,
																Q3,
																number_of_cells_x,
																number_of_cells_y,
																number_of_cells_z,
																cell_overlap)
#	mean_number_of_particles_per_cell::Float64 = 0.0
#	for i = 1:length(cell_lists)
#		mean_number_of_particles_per_cell += length(cell_lists[i])
#	end
#	mean_number_of_particles_per_cell /= length(cell_lists)
#	#	println(cell_lists)
#	println(mean_number_of_particles_per_cell)
#	return

	# Simulate diffusion.
	number_of_workers::Int64 = nworkers() # This is determined by the the '-p' input flag to Julia.
	number_of_diffusers_per_worker::Array{Int64, 1} = convert(Array{Int64, 1}, floor(number_of_diffusers / number_of_workers) * ones(number_of_workers))
	number_of_diffusers_remaining::Int64 = number_of_diffusers - sum(number_of_diffusers_per_worker)
	number_of_diffusers_per_worker[1:number_of_diffusers_remaining] += 1
	println(number_of_diffusers_per_worker)

	output::Array{Float64, 1} = zeros(4*number_of_time_points_coarse+1)
	output = @parallel (+) for current_worker = 1:number_of_workers
		diffuse(
			particle_type,
			R,
			Lx,
			Ly,
			Lz,
			X,
			Y,
			Z,
			Q0,
			Q1,
			Q2,
			Q3,
			A11,
			A12,
			A13,
			A21,
			A22,
			A23,
			A31,
			A32,
			A33,
			D0,
			deltat_coarse,
			number_of_time_points_coarse,
			number_of_time_points_fine_per_coarse,
			number_of_diffusers_per_worker[current_worker],
			boundary_condition,
			cell_lists)
	end
	msd::Array{Float64, 1} = output[1:number_of_time_points_coarse] ./ convert(Float64, 3 * number_of_diffusers)
	msd_x::Array{Float64, 1} = output[number_of_time_points_coarse+1:2*number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	msd_y::Array{Float64, 1} = output[2*number_of_time_points_coarse+1:3*number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	msd_z::Array{Float64, 1} = output[3*number_of_time_points_coarse+1:4*number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	D0_empirical::Float64 = output[end] / (3.0 * convert(Float64, number_of_diffusers * (number_of_time_points_coarse-1) * number_of_time_points_fine_per_coarse) * 2.0 * deltat_fine)
	t_finish_ns::Int64 = convert(Int64, time_ns())
	t_exec::Float64 = convert(Float64, t_finish_ns - t_start_ns) / 1e9

	# Write output.
	diagnostic_diffusion_coefficient_ratio::Float64 = D0_empirical / D0
	write_output(
		output_file_path,
		diagnostic_diffusion_coefficient_ratio,
		deltat_coarse,
		number_of_time_points_coarse,
		msd,
		msd_x,
		msd_y,
		msd_z,
		t_exec)
	println(join(("Output written to ", output_file_path, ".")))
	println("Finished.")

	nothing
end

run_diffusion()
