function write_input(	file_path::String,
						output_generation_path::String,
						diffusion_coefficient::Float64,
						deltat_coarse::Float64,
						number_of_time_points_coarse::Int64,
						number_of_time_points_fine_per_coarse::Int64,
						number_of_diffusers::Int64,
						number_of_cells_x::Int64,
						number_of_cells_y::Int64,
						number_of_cells_z::Int64,
						boundary_condition::String,
						output_file_path::String)

	file_stream::IOStream = open(file_path, "w")

	@printf(file_stream, "%s", "<input_diffusion>\n")

	write_key(file_stream, "output_generation_path", output_generation_path)
	write_key(file_stream, "diffusion_coefficient", diffusion_coefficient)
	write_key(file_stream, "deltat_coarse", deltat_coarse)
	write_key(file_stream, "number_of_time_points_coarse", number_of_time_points_coarse)
	write_key(file_stream, "number_of_time_points_fine_per_coarse", number_of_time_points_fine_per_coarse)
	write_key(file_stream, "number_of_diffusers", number_of_diffusers)
	write_key(file_stream, "number_of_cells_x", number_of_cells_x)
	write_key(file_stream, "number_of_cells_y", number_of_cells_y)
	write_key(file_stream, "number_of_cells_z", number_of_cells_z)
	write_key(file_stream, "boundary_condition", boundary_condition)
	write_key(file_stream, "output_file_path", output_file_path)

	@printf(file_stream, "%s", "</input_diffusion>")

	close(file_stream)

	nothing
end
