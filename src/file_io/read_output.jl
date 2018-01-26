function read_output(file_path::String)
	file_stream::IOStream = open(file_path, "r")
	file_string::String = readstring(file_stream)
	close(file_stream)

	D_empirical::Float64 = read_key(file_string, "D_empirical", Float64)
	t::Array{Float64, 1} = read_key(file_string, "t", Array{Float64, 1})
	msd::Array{Float64, 1} = read_key(file_string, "msd", Array{Float64, 1})
	msd_x::Array{Float64, 1} = read_key(file_string, "msd_x", Array{Float64, 1})
	msd_y::Array{Float64, 1} = read_key(file_string, "msd_y", Array{Float64, 1})
	msd_z::Array{Float64, 1} = read_key(file_string, "msd_z", Array{Float64, 1})
	t_exec::Float64 = read_key(file_string, "t_exec", Float64)

	return (
		D_empirical,
		t,
		msd,
		msd_x,
		msd_y,
		msd_z,
		t_exec)
end
