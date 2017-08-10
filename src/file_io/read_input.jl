function read_input(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	output_generation_path::String = read_key(file_string, "output_generation_path", String)
	diffusion_coefficient::Float64 = read_key(file_string, "diffusion_coefficient", Float64)
	deltat_coarse::Float64 = read_key(file_string, "deltat_coarse", Float64)
	number_of_time_points_coarse::Int64 = read_key(file_string, "number_of_time_points_coarse", Int64)
	number_of_time_points_fine_per_coarse::Int64 = read_key(file_string, "number_of_time_points_fine_per_coarse", Int64)
	number_of_diffusers::Int64 = read_key(file_string, "number_of_diffusers", Int64)
	number_of_cells_x::Int64 = read_key(file_string, "number_of_cells_x", Int64)
	number_of_cells_y::Int64 = read_key(file_string, "number_of_cells_y", Int64)
	number_of_cells_z::Int64 = read_key(file_string, "number_of_cells_z", Int64)
	output_file_path::String = read_key(file_string, "output_file_path", String)

	return (
		output_generation_path,
		diffusion_coefficient,
		deltat_coarse,
		number_of_time_points_coarse,
		number_of_time_points_fine_per_coarse,
		number_of_diffusers,
		number_of_cells_x,
		number_of_cells_y,
		number_of_cells_z,
		output_file_path)
end
