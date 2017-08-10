
include("file_io/read_key.jl")
include("file_io/write_key.jl")

include("file_io/read_input.jl")
include("file_io/read_output_generation.jl")
include("file_io/write_output.jl")

foo = @__FILE__
@eval @everywhere f = $foo
@everywhere println(f)
@everywhere (program_file_dir, program_file_name) = splitdir(f)
@everywhere include(joinpath(program_file_dir, "text_io/print_progress.jl"))
#@everywhere include(joinpath(program_file_dir, "ellipticaldisk/diffuse.jl"))
@everywhere include(joinpath(program_file_dir, "ellipticaldisk/diffuse_parallel.jl"))

function wfrun_diffusion_parallel()
	# Inititalization of random number generation device.
	const random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)

	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)

	# Process input arguments.
	silent_mode::Bool = false
	parallel_mode::Bool = false

	input_file_path::String = ""
	input_file_specified::Bool = false

	number_of_arguments::Int64 = length(ARGS)

	for current_argument = 1:number_of_arguments
		if ARGS[current_argument] == "-silent"
			silent_mode = true
		elseif isfile(ARGS[current_argument])
			input_file_path = ARGS[current_argument]
			input_file_specified = true
		end
	end

	# Print text header.
	if !silent_mode
		print_header()
	end

	# Error checking for input file specification.
	if !input_file_specified
		println("No input file specified or specified input does not exist. Aborting.")
		return nothing
	end

	# Process input file.
	if !silent_mode
		println(join(("Reading input from file ", input_file_path, "...")))
	end
	(output_generation_path::String, D0::Float64, deltat_coarse::Float64, number_of_time_points_coarse::Int64, number_of_time_points_fine_per_coarse::Int64, number_of_diffusers::Int64, number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64, output_file_path::String) = read_input(input_file_path)

	# Process particle system file.
	if !silent_mode
		println(join(("Reading particle system configuration from file ", output_generation_path, "...")))
	end
	(Lx::Float64, Ly::Float64, Lz::Float64, particle_type::String, number_of_particles::Int64, X::Array{Float64, 1}, Y::Array{Float64, 1}, Z::Array{Float64, 1}, THETA1::Array{Float64, 1}, THETA2::Array{Float64, 1}, THETA3::Array{Float64, 1}, R1::Array{Float64, 1}, R2::Array{Float64, 1}) = read_output_generation(output_generation_path)

	# Print simulation stats.
	if !silent_mode
		print_simulation_stats(particle_type, number_of_particles, number_of_diffusers, Lx, Ly, Lz, number_of_cells_x, number_of_cells_y, number_of_cells_z)
	end

	# Simulate diffusion.
	number_of_workers::Int64 = nworkers() # This is determined by the the '-p' input flag to Julia.
	number_of_diffusers_per_worker::Array{Int64, 1} = convert(Array{Int64, 1}, floor(number_of_diffusers / number_of_workers) * ones(number_of_workers))
	number_of_diffusers_remaining::Int64 = number_of_diffusers - sum(number_of_diffusers_per_worker)
	number_of_diffusers_per_worker[1:number_of_diffusers_remaining] += 1
	println(number_of_diffusers_per_worker)

	t_start_ns::Int64 = convert(Int64, time_ns())

	#msd_x = SharedArray(Float64, number_of_time_points_coarse)
	#msd_y = SharedArray(Float64, number_of_time_points_coarse)
	#msd_z = SharedArray(Float64, number_of_time_points_coarse)
	#D0_empirical = SharedArray(Float64, 1)
	#@sync @parallel for current_worker = 1:number_of_workers
		#(msd_x_, msd_y_, msd_z_, D0_empirical_) = diffuse(
			#X,
			#Y,
			#Z,
			#THETA1,
			#THETA2,
			#THETA3,
			#R1,
			#R2,
			#Lx,
			#Ly,
			#Lz,
			#D0,
			#deltat_coarse,
			#number_of_time_points_coarse,
			#number_of_time_points_fine_per_coarse,
			#number_of_diffusers_per_worker[current_worker],
			#number_of_cells_x,
			#number_of_cells_y,
			#number_of_cells_z,
			#silent_mode)
#
		#for i = 1:number_of_time_points_coarse
		#	msd_x[i] = msd_x[i] + convert(Float64, number_of_diffusers_per_worker[current_worker]) / convert(Float64, number_of_diffusers) * msd_x_[i]
		#	msd_y[i] = msd_y[i] + convert(Float64, number_of_diffusers_per_worker[current_worker]) / convert(Float64, number_of_diffusers) * msd_y_[i]
		#	msd_z[i] = msd_z[i] + convert(Float64, number_of_diffusers_per_worker[current_worker]) / convert(Float64, number_of_diffusers) * msd_z_[i]
		#end
		#D0_empirical[1] = D0_empirical[1] + convert(Float64, number_of_diffusers_per_worker[current_worker]) / convert(Float64, number_of_diffusers) * D0_empirical_
	#end
	output::Array{Float64, 1} = zeros(3*number_of_time_points_coarse+1)
	output = @parallel (+) for current_worker = 1:number_of_workers
		diffuse_parallel(
			X,
			Y,
			Z,
			THETA1,
			THETA2,
			THETA3,
			R1,
			R2,
			Lx,
			Ly,
			Lz,
			D0,
			deltat_coarse,
			number_of_time_points_coarse,
			number_of_time_points_fine_per_coarse,
			number_of_diffusers_per_worker[current_worker],
			number_of_cells_x,
			number_of_cells_y,
			number_of_cells_z,
			silent_mode)
	end
	msd_x::Array{Float64, 1} = output[1:number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	msd_y::Array{Float64, 1} = output[number_of_time_points_coarse+1:2*number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	msd_z::Array{Float64, 1} = output[2*number_of_time_points_coarse+1:3*number_of_time_points_coarse] ./ convert(Float64, number_of_diffusers)
	deltat_fine::Float64 = deltat_coarse / convert(Float64, number_of_time_points_fine_per_coarse)
	D0_empirical::Float64 = output[end] / (3.0 * convert(Float64, number_of_diffusers * (number_of_time_points_coarse-1) * number_of_time_points_fine_per_coarse) * 2.0 * deltat_fine)
	t_finish_ns::Int64 = convert(Int64, time_ns())
	t_exec::Float64 = convert(Float64, t_finish_ns - t_start_ns) / 1e9

	# Write output.

	#write_output(output_file_path, D0, convert(Float64, D0_empirical[1]), deltat_coarse, number_of_time_points_coarse, convert(Array{Float64, 1}, msd_x), convert(Array{Float64, 1}, msd_y), convert(Array{Float64, 1}, msd_z), t_exec)
	#write_output(output_file_path, D0, sdata(D0_empirical)[1], deltat_coarse, number_of_time_points_coarse, sdata(msd_x), sdata(msd_y), sdata(msd_z), t_exec)
	write_output(output_file_path, D0, D0_empirical, deltat_coarse, number_of_time_points_coarse, msd_x, msd_y, msd_z, t_exec)

	# Print output information.
	if !silent_mode
		println(join(("Output written to ", output_file_path, ".")))
		println("Finished.")
	end

	nothing
end

wfrun_diffusion_parallel()
