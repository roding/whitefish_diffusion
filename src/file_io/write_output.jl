function write_output(	file_path::String,
						diagnostic_diffusion_coefficient_ratio::Float64,
						deltat_coarse::Float64,
						number_of_time_points_coarse::Int64,
						msd::Array{Float64, 1},
						msd_x::Array{Float64, 1},
						msd_y::Array{Float64, 1},
						msd_z::Array{Float64, 1},
						t_exec::Float64)
	file_stream::IOStream = open(file_path, "w")

	t::Array{Float64, 1} = 0.0:deltat_coarse:(convert(Float64, number_of_time_points_coarse-1) * deltat_coarse)

	@printf(file_stream, "%s", "<output_diffusion>\n")

	write_key(file_stream, "diagnostic_diffusion_coefficient_ratio", diagnostic_diffusion_coefficient_ratio)
	write_key(file_stream, "time", t)
	write_key(file_stream, "mean_square_displacement", msd)
	write_key(file_stream, "mean_square_displacement_x", msd_x)
	write_key(file_stream, "mean_square_displacement_y", msd_y)
	write_key(file_stream, "mean_square_displacement_z", msd_z)
	write_key(file_stream, "execution_time", t_exec)

	@printf(file_stream, "%s", "</output_diffusion>")

	close(file_stream)

	nothing
end
