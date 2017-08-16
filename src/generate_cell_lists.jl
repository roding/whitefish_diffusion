function generate_cell_lists(	particle_type::String,
								R::Array{Float64, 2},
								Lx::Float64,
								Ly::Float64,
								Lz::Float64,
								X::Array{Float64, 1},
								Y::Array{Float64, 1},
								Z::Array{Float64, 1},
								Q0::Array{Float64, 1},
								Q1::Array{Float64, 1},
								Q2::Array{Float64, 1},
								Q3::Array{Float64, 1},
								number_of_cells_x::Int64,
								number_of_cells_y::Int64,
								number_of_cells_z::Int64,
								cell_overlap::Float64)

	# Number of particles.
	number_of_particles::Int64 = length(X)

	# Create cell dimension data structures.
	cell_bounds_x::Array{Float64, 1} = linspace(0.0, Lx, number_of_cells_x + 1)
	lbx_cell::Array{Float64, 1} = cell_bounds_x[1:end-1] - cell_overlap
	ubx_cell::Array{Float64, 1} = cell_bounds_x[2:end] + cell_overlap

	cell_bounds_y::Array{Float64, 1} = linspace(0.0, Ly, number_of_cells_y + 1)
	lby_cell::Array{Float64, 1} = cell_bounds_y[1:end-1] - cell_overlap
	uby_cell::Array{Float64, 1} = cell_bounds_y[2:end] + cell_overlap

	cell_bounds_z::Array{Float64, 1} = linspace(0.0, Lz, number_of_cells_z + 1)
	lbz_cell::Array{Float64, 1} = cell_bounds_z[1:end-1] - cell_overlap
	ubz_cell::Array{Float64, 1} = cell_bounds_z[2:end] + cell_overlap

	# Cell lists data structure.
	cell_lists = Array{Array{Int64, 1}}(number_of_cells_x, number_of_cells_y, number_of_cells_z)
	for current_cell_x = 1:number_of_cells_x
		for current_cell_y = 1:number_of_cells_y
			for current_cell_z = 1:number_of_cells_z
				cell_lists[current_cell_x, current_cell_y, current_cell_z] = Array{Int64}(0)
			end
		end
	end

	# Compute cell lists.
	a11::Float64 = 0.0
	a12::Float64 = 0.0
	a13::Float64 = 0.0
	a21::Float64 = 0.0
	a22::Float64 = 0.0
	a23::Float64 = 0.0
	a31::Float64 = 0.0
	a32::Float64 = 0.0
	a33::Float64 = 0.0
	b11::Float64 = 0.0
	b12::Float64 = 0.0
	b13::Float64 = 0.0
	b21::Float64 = 0.0
	b22::Float64 = 0.0
	b23::Float64 = 0.0
	b31::Float64 = 0.0
	b32::Float64 = 0.0
	b33::Float64 = 0.0
	xAB::Float64 = 0.0
	yAB::Float64 = 0.0
	zAB::Float64 = 0.0

	if particle_type == "sphere"

	elseif particle_type == "ellipse"

	elseif particle_type == "ellipsoid"

	elseif particle_type == "cuboid"
		(b11, b12, b13, b21, b22, b23, b31, b32, b33) = rotation_matrix(1.0, 0.0, 0.0, 0.0) # Rotation matrix of (axis-aligned) cell.
		for current_cell_x = 1:number_of_cells_x
			for current_cell_y = 1:number_of_cells_y
				for current_cell_z = 1:number_of_cells_z
					for current_particle = 1:number_of_particles
						xAB = signed_distance_mod(X[current_particle], 0.5 * (lbx_cell[current_cell_x] + ubx_cell[current_cell_x]), Lx)
						yAB = signed_distance_mod(Y[current_particle], 0.5 * (lby_cell[current_cell_y] + uby_cell[current_cell_y]), Ly)
						zAB = signed_distance_mod(Z[current_particle], 0.5 * (lbz_cell[current_cell_z] + ubz_cell[current_cell_z]), Lz)
						(a11, a12, a13, a21, a22, a23, a31, a32, a33) = rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
						if overlap_cuboid_binary(xAB, yAB, zAB, a11, a12, a13, a21, a22, a23, a31, a32, a33, b11, b12, b13, b21, b22, b23, b31, b32, b33,
							R[current_particle, 1], R[current_particle, 2], R[current_particle, 3],
							0.5 * (ubx_cell[current_cell_x] - lbx_cell[current_cell_x]),
							0.5 * (uby_cell[current_cell_y] - lby_cell[current_cell_y]),
							0.5 * (ubz_cell[current_cell_z] - lbz_cell[current_cell_z])) == 1.0

							push!(cell_lists[current_cell_x, current_cell_y, current_cell_z], current_particle)
						end
					end
				end
			end
		end
	end



#					if intersect_box_box(	lbx_particle[current_particle],
#											ubx_particle[current_particle],
#											lby_particle[current_particle],
#											uby_particle[current_particle],
#											lbz_particle[current_particle],
#											ubz_particle[current_particle],
#											lbx_cell[current_cell_x],
#											ubx_cell[current_cell_x],
#											lby_cell[current_cell_y],
#											uby_cell[current_cell_y],
#											lbz_cell[current_cell_z],
#											ubz_cell[current_cell_z],
#											Lx,
#											Ly,
#											Lz)



	# Compute bounding box for all particles.
#	lbx_particle::Array{Float64, 1} = zeros(number_of_particles)
#	ubx_particle::Array{Float64, 1} = zeros(number_of_particles)
#	lby_particle::Array{Float64, 1} = zeros(number_of_particles)
#	uby_particle::Array{Float64, 1} = zeros(number_of_particles)
#	lbz_particle::Array{Float64, 1} = zeros(number_of_particles)
#	ubz_particle::Array{Float64, 1} = zeros(number_of_particles)
#
#	for current_particle = 1:number_of_particles
#		(	lbx_particle[current_particle],
#			ubx_particle[current_particle],
#			lby_particle[current_particle],
#			uby_particle[current_particle],
#			lbz_particle[current_particle],
#			ubz_particle[current_particle]) = axis_aligned_bounding_box(particle_type,
#																		R[current_particle, :],
#																		X[current_particle],
#																		Y[current_particle],
#																		Z[current_particle],
#																		Q0[current_particle],
#																		Q1[current_particle],
#																		Q2[current_particle],
#																		Q3[current_particle])
#	end





	return cell_lists
end
