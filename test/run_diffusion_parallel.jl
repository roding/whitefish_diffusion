include("../src/file_io/write_xml_key.jl")
include("../src/file_io/write_xml_input.jl")
include("../src/file_io/read_xml_output_generation.jl")

function run_diffusion_parallel(number_of_workers::Int64)
	# Inititalization of random number generation device.
	random_seed::Int64 = convert(Int64, time_ns())
	srand(random_seed)

	# Change current folder to the folder where this script lies.
	(program_file_dir::String, program_file_name::String) = splitdir(PROGRAM_FILE)
	program_file_dir = abspath(program_file_dir)
	cd(program_file_dir)

	# Output directory.
	output_dir::String = abspath(joinpath(program_file_dir, "../../output"))
	if !isdir(output_dir)
        mkdir(output_dir)
    end

	# Diffusion parameters.
	D0::Float64 = 1.0
	deltat_coarse::Float64 = 0.25
	number_of_time_points_coarse::Int64 = 50000
	number_of_time_points_fine_per_coarse::Int64 = 100
	number_of_diffusers::Int64 = 100000#100000
	number_of_cells_x::Int64 = 10
	number_of_cells_y::Int64 = 10
	number_of_cells_z::Int64 = 10
	output_file_path::String = joinpath(output_dir, "output_diffusion.xml")

	# Write input file for diffusion.
	output_generation_path::String = joinpath(output_dir, "output_generation.xml")
	input_file_path::String = joinpath(output_dir, "input_diffusion.xml")
	write_xml_input(input_file_path, output_generation_path, D0, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z, output_file_path)

	# Run diffusion.
	program_path::String = abspath("../src/wfrun_diffusion_parallel.jl")
	cmd::Cmd = `julia -p $number_of_workers $program_path $input_file_path`
	run(cmd)

	# Exit.
	nothing

end

run_diffusion_parallel(88)
