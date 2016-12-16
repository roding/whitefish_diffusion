function read_xml_test()
	file_name::String = "particle_system.xml"
	file_stream::IOStream = open(file_name, "r")
	#file_lines::Array{String, 1} = readlines(file_stream)
	#println(file_lines)
	file_string::String = readstring(file_stream)
	close(file_stream)
	
	# Read simulation domain size in x direction.
	ind_before = search(file_string, "<domain_size_x>")
	ind_after = search(file_string, "</domain_size_x>")
	ind = ind_before[end]+1:ind_after[1]-1
	Lx::Float64 = parse(file_string[ind])
	
	# Read simulation domain size in y direction.
	ind_before = search(file_string, "<domain_size_y>")
	ind_after = search(file_string, "</domain_size_y>")
	ind = ind_before[end]+1:ind_after[1]-1
	Ly::Float64 = parse(file_string[ind])
	
	# Read simulation domain size in z direction.
	ind_before = search(file_string, "<domain_size_z>")
	ind_after = search(file_string, "</domain_size_z>")
	ind = ind_before[end]+1:ind_after[1]-1
	Lz::Float64 = parse(file_string[ind])
	
	# Read particle type.
	ind_before = search(file_string, "<type>")
	ind_after = search(file_string, "</type>")
	ind = ind_before[end]+1:ind_after[1]-1
	particle_category::Symbol = Symbol(file_string[ind])
	
	# Read number of particles.
	ind_before = search(file_string, "<number_of_particles>")
	ind_after = search(file_string, "</number_of_particles>")
	ind = ind_before[end]+1:ind_after[1]-1
	number_of_particles::Int64 = parse(file_string[ind])
	
	# For now we assume particles are ellipical disks.
	X = Array(Float64, number_of_particles)
	Y = Array(Float64, number_of_particles)
	Z = Array(Float64, number_of_particles)
	THETA1 = Array(Float64, number_of_particles)
	THETA2 = Array(Float64, number_of_particles)
	THETA3 = Array(Float64, number_of_particles)
	R1 = Array(Float64, number_of_particles)
	R2 = Array(Float64, number_of_particles)
	ind_before = search(file_string, "<parameters>")
	ind_after = search(file_string, "</parameters>")
	ind = ind_before[end]+1:ind_after[1]-1
	parameters_string_array::Array{String, 1} = split(file_string[ind], ",")
	for current_particle = 1:number_of_particles
		
	
	parameters_array = Array(Float64, number_of_particles * 8) = 
	
	
	
	
	
	
	println(Lx)
	#stream_ = IOBuffer(string)
	
	
	
	
	
	println(ind_before)
	println(ind_after)





	nothing
end

read_xml_test()