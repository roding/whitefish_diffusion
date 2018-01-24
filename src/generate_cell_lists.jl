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
	xAB::Float64 = 0.0
	yAB::Float64 = 0.0
	zAB::Float64 = 0.0

	if particle_type == "sphere"
		# For now we do a simple AABB-AABB intersection test.
		# Should implement the exact test later because there will be some false positives now.
		for current_cell_x = 1:number_of_cells_x
			for current_cell_y = 1:number_of_cells_y
				for current_cell_z = 1:number_of_cells_z
					for current_particle = 1:number_of_particles
						xAB = signed_distance_mod(X[current_particle], 0.5 * (lbx_cell[current_cell_x] + ubx_cell[current_cell_x]), Lx)
						yAB = signed_distance_mod(Y[current_particle], 0.5 * (lby_cell[current_cell_y] + uby_cell[current_cell_y]), Ly)
						zAB = signed_distance_mod(Z[current_particle], 0.5 * (lbz_cell[current_cell_z] + ubz_cell[current_cell_z]), Lz)
						if overlap_cuboid_binary(xAB, yAB, zAB, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
							R[current_particle, 1], R[current_particle, 1], R[current_particle, 1], # SIC!
							0.5 * (ubx_cell[current_cell_x] - lbx_cell[current_cell_x]),
							0.5 * (uby_cell[current_cell_y] - lby_cell[current_cell_y]),
							0.5 * (ubz_cell[current_cell_z] - lbz_cell[current_cell_z])) == 1.0

							push!(cell_lists[current_cell_x, current_cell_y, current_cell_z], current_particle)
						end
					end
				end
			end
		end
	elseif particle_type == "ellipse"
		# For now we approximate the ellipse with its OBB assuming a small finite thickness.
		# Should implement the exact test later because there will be some false positives now.
		for current_cell_x = 1:number_of_cells_x
			for current_cell_y = 1:number_of_cells_y
				for current_cell_z = 1:number_of_cells_z
					for current_particle = 1:number_of_particles
						xAB = signed_distance_mod(X[current_particle], 0.5 * (lbx_cell[current_cell_x] + ubx_cell[current_cell_x]), Lx)
						yAB = signed_distance_mod(Y[current_particle], 0.5 * (lby_cell[current_cell_y] + uby_cell[current_cell_y]), Ly)
						zAB = signed_distance_mod(Z[current_particle], 0.5 * (lbz_cell[current_cell_z] + ubz_cell[current_cell_z]), Lz)
						(a11, a12, a13, a21, a22, a23, a31, a32, a33) = rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
						if overlap_cuboid_binary(xAB, yAB, zAB, a11, a12, a13, a21, a22, a23, a31, a32, a33, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
							R[current_particle, 1], R[current_particle, 2], 0.1, # 0.1 is the artificial finite thickness of the OBB.
							0.5 * (ubx_cell[current_cell_x] - lbx_cell[current_cell_x]),
							0.5 * (uby_cell[current_cell_y] - lby_cell[current_cell_y]),
							0.5 * (ubz_cell[current_cell_z] - lbz_cell[current_cell_z])) == 1.0

							push!(cell_lists[current_cell_x, current_cell_y, current_cell_z], current_particle)
						end
					end
				end
			end
		end
	elseif particle_type == "ellipsoid"
		# For now we approximate the ellipsoid with its OBB i.e. the bounding cuboid with same semi-axes.
		# Should implement the exact test later because there will be some false positives now.
		for current_cell_x = 1:number_of_cells_x
			for current_cell_y = 1:number_of_cells_y
				for current_cell_z = 1:number_of_cells_z
					for current_particle = 1:number_of_particles
						xAB = signed_distance_mod(X[current_particle], 0.5 * (lbx_cell[current_cell_x] + ubx_cell[current_cell_x]), Lx)
						yAB = signed_distance_mod(Y[current_particle], 0.5 * (lby_cell[current_cell_y] + uby_cell[current_cell_y]), Ly)
						zAB = signed_distance_mod(Z[current_particle], 0.5 * (lbz_cell[current_cell_z] + ubz_cell[current_cell_z]), Lz)
						(a11, a12, a13, a21, a22, a23, a31, a32, a33) = rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
						if overlap_cuboid_binary(xAB, yAB, zAB, a11, a12, a13, a21, a22, a23, a31, a32, a33, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
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
	elseif particle_type == "cuboid"
		for current_cell_x = 1:number_of_cells_x
			for current_cell_y = 1:number_of_cells_y
				for current_cell_z = 1:number_of_cells_z
					for current_particle = 1:number_of_particles
						xAB = signed_distance_mod(X[current_particle], 0.5 * (lbx_cell[current_cell_x] + ubx_cell[current_cell_x]), Lx)
						yAB = signed_distance_mod(Y[current_particle], 0.5 * (lby_cell[current_cell_y] + uby_cell[current_cell_y]), Ly)
						zAB = signed_distance_mod(Z[current_particle], 0.5 * (lbz_cell[current_cell_z] + ubz_cell[current_cell_z]), Lz)
						(a11, a12, a13, a21, a22, a23, a31, a32, a33) = rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
						if overlap_cuboid_binary(xAB, yAB, zAB, a11, a12, a13, a21, a22, a23, a31, a32, a33, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
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
	elseif particle_type == "superellipsoid"
		# For now we approximate the ellipsoid with its OBB i.e. the bounding cuboid with the same semi-axes.
		# Should implement the exact test later because there will be some false positives now.
		for current_cell_x = 1:number_of_cells_x
			for current_cell_y = 1:number_of_cells_y
				for current_cell_z = 1:number_of_cells_z
					for current_particle = 1:number_of_particles
						xAB = signed_distance_mod(X[current_particle], 0.5 * (lbx_cell[current_cell_x] + ubx_cell[current_cell_x]), Lx)
						yAB = signed_distance_mod(Y[current_particle], 0.5 * (lby_cell[current_cell_y] + uby_cell[current_cell_y]), Ly)
						zAB = signed_distance_mod(Z[current_particle], 0.5 * (lbz_cell[current_cell_z] + ubz_cell[current_cell_z]), Lz)
						(a11, a12, a13, a21, a22, a23, a31, a32, a33) = rotation_matrix(Q0[current_particle], Q1[current_particle], Q2[current_particle], Q3[current_particle])
						if overlap_cuboid_binary(xAB, yAB, zAB, a11, a12, a13, a21, a22, a23, a31, a32, a33, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
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

	return cell_lists
end
