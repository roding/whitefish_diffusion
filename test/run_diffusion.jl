include("../src/io/write_xml_key.jl")
include("../src/io/write_xml_input.jl")
include("../src/io/write_xml_input_generation.jl")

function run_diffusion()
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
	deltat_coarse::Float64 = 5.0
	number_of_time_points_coarse::Int64 = 20000
	number_of_time_points_fine_per_coarse::Int64 = 5#125
	number_of_diffusers::Int64 = 256000
	number_of_cells_x::Int64 = 10
	number_of_cells_y::Int64 = 10
	number_of_cells_z::Int64 = 10
	output_file_path::String = joinpath(output_dir, "output_diffusion.xml")
	
	# Write input file for diffusion.
	output_generation_path::String = joinpath(output_dir, "output_generation.xml")
	input_file_path::String = joinpath(output_dir, "input_diffusion.xml")
	write_xml_input_diffusion(input_file_path, output_generation_path, D0, deltat_coarse, number_of_time_points_coarse, number_of_time_points_fine_per_coarse, number_of_diffusers, number_of_cells_x, number_of_cells_y, number_of_cells_z, output_file_path)
	
	# Run diffusion.
	program_path::String = abspath("../src/wfrun_diffusion.jl")
	cmd::Cmd = `julia $program_path $input_file_path`
	run(cmd_diffusion)
	
	# Exit.
	nothing
	
end

run_diffusion()