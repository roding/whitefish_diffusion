function diffuse(	particle_type_ss::String,
					R_ss::Array{Float64, 2},
					X_ss::Array{Float64, 1},
					Y_ss::Array{Float64, 1},
					Z_ss::Array{Float64, 1},
					Q0_ss::Array{Float64, 1},
					Q1_ss::Array{Float64, 1},
					Q2_ss::Array{Float64, 1},
					Q3_ss::Array{Float64, 1},
					A11_ss::Array{Float64, 1},
					A12_ss::Array{Float64, 1},
					A13_ss::Array{Float64, 1},
					A21_ss::Array{Float64, 1},
					A22_ss::Array{Float64, 1},
					A23_ss::Array{Float64, 1},
					A31_ss::Array{Float64, 1},
					A32_ss::Array{Float64, 1},
					A33_ss::Array{Float64, 1},
					cell_lists_ss::Array{Array{Int64, 1}, 3},
					particle_type_dm::String,
					R_dm::Array{Float64, 2},
					X_dm::Array{Float64, 1},
					Y_dm::Array{Float64, 1},
					Z_dm::Array{Float64, 1},
					Q0_dm::Array{Float64, 1},
					Q1_dm::Array{Float64, 1},
					Q2_dm::Array{Float64, 1},
					Q3_dm::Array{Float64, 1},
					A11_dm::Array{Float64, 1},
					A12_dm::Array{Float64, 1},
					A13_dm::Array{Float64, 1},
					A21_dm::Array{Float64, 1},
					A22_dm::Array{Float64, 1},
					A23_dm::Array{Float64, 1},
					A31_dm::Array{Float64, 1},
					A32_dm::Array{Float64, 1},
					A33_dm::Array{Float64, 1},
					cell_lists_dm::Array{Array{Int64, 1}, 3},
					Lx::Float64,
					Ly::Float64,
					Lz::Float64,
					D0::Float64,
					D1::Float64,
					deltat_coarse::Float64,
					number_of_time_points_coarse::Int64,
					number_of_time_points_fine_per_coarse::Int64,
					number_of_diffusers::Int64,
					boundary_condition::String)

	# Standard deviations of Gaussian jumps.
	deltat_fine::Float64 = deltat_coarse / convert(Float64, number_of_time_points_fine_per_coarse)
	sigma0::Float64 = sqrt(2.0 * D0 * deltat_fine)
	sigma1::Float64 = sqrt(2.0 * D1 * deltat_fine)
	sigma::Float64 = 0.0

	# Number of cells (number of cells must be the same for both cell lists)
	(number_of_cells_x::Int64, number_of_cells_y::Int64, number_of_cells_z::Int64) = size(cell_lists_ss)

	current_particle_in_cell::Int64 = 0
	is_proposed_position_ok::Bool = true
	is_in_phase0::Bool = true

	x::Float64 = 0.0
	y::Float64 = 0.0
	z::Float64 = 0.0
	x_abs::Float64 = 0.0
	y_abs::Float64 = 0.0
	z_abs::Float64 = 0.0
	x_star::Float64 = 0.0
	y_star::Float64 = 0.0
	z_star::Float64 = 0.0
	deltax::Float64 = 0.0
	deltay::Float64 = 0.0
	deltaz::Float64 = 0.0
	vx::Float64 = 0.0
	vy::Float64 = 0.0
	vz::Float64 = 0.0
	vx_star::Float64 = 0.0
	vy_star::Float64 = 0.0
	vz_star::Float64 = 0.0
	w1::Float64 = 0.0
	w1_star::Float64 = 0.0
	w2::Float64 = 0.0
	w2_star::Float64 = 0.0
	w3::Float64 = 0.0
	w3_star::Float64 = 0.0
	w1_intersection::Float64 = 0.0
	w2_intersection::Float64 = 0.0
	alpha::Float64 = 0.0
	current_cell_x::Int64 = 0
	current_cell_y::Int64 = 0
	current_cell_z::Int64 = 0
	number_of_particles_current_cell::Int64 = 0
	trajectory_x::Array{Float64} = zeros(number_of_time_points_coarse)
	trajectory_y::Array{Float64} = zeros(number_of_time_points_coarse)
	trajectory_z::Array{Float64} = zeros(number_of_time_points_coarse)
	ssd::Array{Float64} = zeros(number_of_time_points_coarse)
	ssd_x::Array{Float64} = zeros(number_of_time_points_coarse)
	ssd_y::Array{Float64} = zeros(number_of_time_points_coarse)
	ssd_z::Array{Float64} = zeros(number_of_time_points_coarse)
	D_empirical::Float64 = 0.0

	t_start_diffusion::Float64 = convert(Float64, time_ns()) / 1e9
	t_elapsed_diffusion::Float64 = 0.0

	chunk::Int64 = 0
	for current_diffuser = 1:number_of_diffusers
		if mod(current_diffuser, 10) == 0
			println(current_diffuser)
		end

		# Pick random, valid initial position.
		if particle_type_ss == "sphere"
			is_proposed_position_ok = false
			while !is_proposed_position_ok
				x = Lx * rand()
				y = Ly * rand()
				z = Lz * rand()

				current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
				current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
				current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
				number_of_particles_current_cell = length(cell_lists_ss[current_cell_x, current_cell_y, current_cell_z])
				current_cell_list_ss = cell_lists_ss[current_cell_x, current_cell_y, current_cell_z]

				is_proposed_position_ok = true
				current_particle_in_cell = 0
				while current_particle_in_cell < number_of_particles_current_cell && is_proposed_position_ok
					current_particle_in_cell += 1
					vx = signed_distance_mod(x, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
					vy = signed_distance_mod(y, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
					vz = signed_distance_mod(z, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
					if vx^2 + vy^2 + vz^2 <= R_ss[current_cell_list_ss[current_particle_in_cell], 1]^2
						is_proposed_position_ok = false
					end
				end
			end
		elseif particle_type_ss == "ellipse"
			# By definition, any point in the simulation domain is outside of the ellipses w.p. 1.
			x = Lx * rand()
			y = Ly * rand()
			z = Lz * rand()
		elseif particle_type_ss == "ellipsoid"
			is_proposed_position_ok = false
			while !is_proposed_position_ok
				x = Lx * rand()
				y = Ly * rand()
				z = Lz * rand()

				current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
				current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
				current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
				number_of_particles_current_cell = length(cell_lists_ss[current_cell_x, current_cell_y, current_cell_z])
				current_cell_list_ss = cell_lists_ss[current_cell_x, current_cell_y, current_cell_z]

				is_proposed_position_ok = true
				current_particle_in_cell = 0
				while current_particle_in_cell < number_of_particles_current_cell && is_proposed_position_ok
					current_particle_in_cell += 1
					vx = signed_distance_mod(x, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
					vy = signed_distance_mod(y, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
					vz = signed_distance_mod(z, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)

					if vx * (A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz) + vy * (A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz) + vz * (A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz) <= 1.0
						is_proposed_position_ok = false
					end
				end
			end
		elseif particle_type_ss == "cuboid"
			is_proposed_position_ok = false
			while !is_proposed_position_ok
				x = Lx * rand()
				y = Ly * rand()
				z = Lz * rand()

				current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
				current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
				current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
				number_of_particles_current_cell = length(cell_lists_ss[current_cell_x, current_cell_y, current_cell_z])
				current_cell_list_ss = cell_lists_ss[current_cell_x, current_cell_y, current_cell_z]

				is_proposed_position_ok = true
				current_particle_in_cell = 0
				while current_particle_in_cell < number_of_particles_current_cell && is_proposed_position_ok
					current_particle_in_cell += 1
					vx = signed_distance_mod(x, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
					vy = signed_distance_mod(y, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
					vz = signed_distance_mod(z, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
					(vx, vy, vz) = (A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz,
									A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz,
									A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz)

					if abs(vx) <= R_ss[current_cell_list_ss[current_particle_in_cell], 1] && abs(vy) <= R_ss[current_cell_list_ss[current_particle_in_cell], 2] && abs(vz) <= R_ss[current_cell_list_ss[current_particle_in_cell], 3]
						is_proposed_position_ok = false
					end
				end
			end
		elseif particle_type_ss == "superellipsoid"
			is_proposed_position_ok = false
			while !is_proposed_position_ok
				x = Lx * rand()
				y = Ly * rand()
				z = Lz * rand()

				current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
				current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
				current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
				number_of_particles_current_cell = length(cell_lists_ss[current_cell_x, current_cell_y, current_cell_z])
				current_cell_list_ss = cell_lists_ss[current_cell_x, current_cell_y, current_cell_z]

				is_proposed_position_ok = true
				current_particle_in_cell = 0
				while current_particle_in_cell < number_of_particles_current_cell && is_proposed_position_ok
					current_particle_in_cell += 1
					vx = signed_distance_mod(x, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
					vy = signed_distance_mod(y, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
					vz = signed_distance_mod(z, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
					(vx, vy, vz) = (A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz,
									A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz,
									A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz)

					if (abs(vx)/R_ss[current_cell_list_ss[current_particle_in_cell], 1])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] +
						(abs(vy)/R_ss[current_cell_list_ss[current_particle_in_cell], 2])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] +
						(abs(vz)/R_ss[current_cell_list_ss[current_particle_in_cell], 3])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] <= 1.0
						is_proposed_position_ok = false
					end
				end
			end
		end

		x_abs = x
		y_abs = y
		z_abs = z

		trajectory_x[1] = x_abs
		trajectory_y[1] = y_abs
		trajectory_z[1] = z_abs

		# Starting diffusion.
		for current_time_coarse = 2:number_of_time_points_coarse
			for current_time_fine = 1:number_of_time_points_fine_per_coarse
				# Diffusion map.
				if particle_type_dm == "sphere"
					current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
					current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
					current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
					number_of_particles_current_cell = length(cell_lists_dm[current_cell_x, current_cell_y, current_cell_z])
					current_cell_list_dm = cell_lists_dm[current_cell_x, current_cell_y, current_cell_z]

					is_in_phase0 = true
					current_particle_in_cell = 0
					while current_particle_in_cell < number_of_particles_current_cell && is_in_phase0
						current_particle_in_cell += 1
						vx = signed_distance_mod(x, X_dm[current_cell_list_dm[current_particle_in_cell]], Lx)
						vy = signed_distance_mod(y, Y_dm[current_cell_list_dm[current_particle_in_cell]], Ly)
						vz = signed_distance_mod(z, Z_dm[current_cell_list_dm[current_particle_in_cell]], Lz)
						if vx^2 + vy^2 + vz^2 <= R_dm[current_cell_list_dm[current_particle_in_cell], 1]^2
							is_in_phase0 = false
						end
					end
				elseif particle_type_dm == "ellipse"
					# By definition, any point in the simulation domain is in phase0. Ellipses should not be even be used here of course.
					is_in_phase0 = true
				elseif particle_type_dm == "ellipsoid"
					current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
					current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
					current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
					number_of_particles_current_cell = length(cell_lists_dm[current_cell_x, current_cell_y, current_cell_z])
					current_cell_list_dm = cell_lists_dm[current_cell_x, current_cell_y, current_cell_z]

					is_in_phase0 = true
					current_particle_in_cell = 0
					while current_particle_in_cell < number_of_particles_current_cell && is_in_phase0
						current_particle_in_cell += 1
						vx = signed_distance_mod(x, X_dm[current_cell_list_dm[current_particle_in_cell]], Lx)
						vy = signed_distance_mod(y, Y_dm[current_cell_list_dm[current_particle_in_cell]], Ly)
						vz = signed_distance_mod(z, Z_dm[current_cell_list_dm[current_particle_in_cell]], Lz)

						if vx * (A11_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A12_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A13_dm[current_cell_list_dm[current_particle_in_cell]] * vz) + vy * (A21_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A22_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A23_dm[current_cell_list_dm[current_particle_in_cell]] * vz) + vz * (A31_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A32_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A33_dm[current_cell_list_dm[current_particle_in_cell]] * vz) <= 1.0
							is_in_phase0 = false
						end
					end
				elseif particle_type_dm == "cuboid"
					current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
					current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
					current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
					number_of_particles_current_cell = length(cell_lists_dm[current_cell_x, current_cell_y, current_cell_z])
					current_cell_list_dm = cell_lists_dm[current_cell_x, current_cell_y, current_cell_z]

					is_in_phase0 = true
					current_particle_in_cell = 0
					while current_particle_in_cell < number_of_particles_current_cell && is_in_phase0
						current_particle_in_cell += 1
						vx = signed_distance_mod(x, X_dm[current_cell_list_dm[current_particle_in_cell]], Lx)
						vy = signed_distance_mod(y, Y_dm[current_cell_list_dm[current_particle_in_cell]], Ly)
						vz = signed_distance_mod(z, Z_dm[current_cell_list_dm[current_particle_in_cell]], Lz)
						(vx, vy, vz) = (A11_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A12_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A13_dm[current_cell_list_dm[current_particle_in_cell]] * vz,
										A21_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A22_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A23_dm[current_cell_list_dm[current_particle_in_cell]] * vz,
										A31_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A32_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A33_dm[current_cell_list_dm[current_particle_in_cell]] * vz)

						if abs(vx) <= R_dm[current_cell_list_dm[current_particle_in_cell], 1] && abs(vy) <= R_dm[current_cell_list_dm[current_particle_in_cell], 2] && abs(vz) <= R_dm[current_cell_list_dm[current_particle_in_cell], 3]
							is_in_phase0 = false
						end
					end
				elseif particle_type_dm == "superellipsoid"
					current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
					current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
					current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
					number_of_particles_current_cell = length(cell_lists_dm[current_cell_x, current_cell_y, current_cell_z])
					current_cell_list_dm = cell_lists_dm[current_cell_x, current_cell_y, current_cell_z]

					is_in_phase0 = true
					current_particle_in_cell = 0
					while current_particle_in_cell < number_of_particles_current_cell && is_in_phase0
						current_particle_in_cell += 1
						vx = signed_distance_mod(x, X_dm[current_cell_list_dm[current_particle_in_cell]], Lx)
						vy = signed_distance_mod(y, Y_dm[current_cell_list_dm[current_particle_in_cell]], Ly)
						vz = signed_distance_mod(z, Z_dm[current_cell_list_dm[current_particle_in_cell]], Lz)
						(vx, vy, vz) = (A11_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A12_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A13_dm[current_cell_list_dm[current_particle_in_cell]] * vz,
										A21_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A22_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A23_dm[current_cell_list_dm[current_particle_in_cell]] * vz,
										A31_dm[current_cell_list_dm[current_particle_in_cell]] * vx + A32_dm[current_cell_list_dm[current_particle_in_cell]] * vy + A33_dm[current_cell_list_dm[current_particle_in_cell]] * vz)

						if (abs(vx)/R_dm[current_cell_list_dm[current_particle_in_cell], 1])^R_dm[current_cell_list_dm[current_particle_in_cell], 4] +
							(abs(vy)/R_dm[current_cell_list_dm[current_particle_in_cell], 2])^R_dm[current_cell_list_dm[current_particle_in_cell], 4] +
							(abs(vz)/R_dm[current_cell_list_dm[current_particle_in_cell], 3])^R_dm[current_cell_list_dm[current_particle_in_cell], 4] <= 1.0
							is_in_phase0 = false
						end
					end
				end
				if is_in_phase0
					sigma = sigma0
				else
					sigma = sigma1
				end

				# Solid structure.
				current_cell_x = convert(Int64, ceil(x / Lx * convert(Float64, number_of_cells_x)))
				current_cell_y = convert(Int64, ceil(y / Ly * convert(Float64, number_of_cells_y)))
				current_cell_z = convert(Int64, ceil(z / Lz * convert(Float64, number_of_cells_z)))
				number_of_particles_current_cell = length(cell_lists_ss[current_cell_x, current_cell_y, current_cell_z])
				current_cell_list_ss = cell_lists_ss[current_cell_x, current_cell_y, current_cell_z]
				if boundary_condition == "single-rejection"
					deltax = sigma * randn()
					deltay = sigma * randn()
					deltaz = sigma * randn()

					x_star = position_mod(x + deltax, Lx)
					y_star = position_mod(y + deltay, Ly)
					z_star = position_mod(z + deltaz, Lz)

					# Check for diffuser-particle intersections.
					current_particle_in_cell = 0
					is_proposed_position_ok = true
					while current_particle_in_cell < number_of_particles_current_cell && is_proposed_position_ok
						current_particle_in_cell += 1

						if particle_type_ss == "sphere"
							vx_star = signed_distance_mod(x_star, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
							vy_star = signed_distance_mod(y_star, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
							vz_star = signed_distance_mod(z_star, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)

							if vx_star^2 + vy_star^2 + vz_star^2 <= R_ss[current_cell_list_ss[current_particle_in_cell], 1]^2
								is_proposed_position_ok = false
							end
						elseif particle_type_ss == "ellipse"
							# Make sure we can detect intersection by ensuring that v and v_star are close, not just close mod L.
							# Ellipse is a bit special in this regard.

							# Observe that the arguments to signed_distance_mod are shifted. Calling them
							# the same way as in other place will lead to wrong sign of distances.
							vx = signed_distance_mod(X_ss[current_cell_list_ss[current_particle_in_cell]], x, Lx)
							vy = signed_distance_mod(Y_ss[current_cell_list_ss[current_particle_in_cell]], y, Ly)
							vz = signed_distance_mod(Z_ss[current_cell_list_ss[current_particle_in_cell]], z, Lz)
							vx_star = vx + deltax
							vy_star = vy + deltay
							vz_star = vz + deltaz
							(vx, vy, vz) = (				A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx +	A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz,
															A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz,
															A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz)
							(vx_star, vy_star, vz_star) = (	A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
															A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
															A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star)

							alpha = - vz_star / (vz - vz_star)
							if 0.0 < alpha < 1.0
								w1_intersection = alpha * vx + (1.0 - alpha) * vx_star
								w2_intersection = alpha * vy + (1.0 - alpha) * vy_star
								if (w1_intersection/R_ss[current_cell_list_ss[current_particle_in_cell], 1])^2 + (w2_intersection/R_ss[current_cell_list_ss[current_particle_in_cell], 2])^2 < 1.0
									is_proposed_position_ok = false
								end
							end
						elseif particle_type_ss == "ellipsoid"
							vx_star = signed_distance_mod(x_star, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
							vy_star = signed_distance_mod(y_star, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
							vz_star = signed_distance_mod(z_star, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
							if vx_star * (A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star) + vy_star * (A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star) + vz_star * (A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star) <= 1.0
								is_proposed_position_ok = false
							end
						elseif particle_type_ss == "cuboid"
							vx_star = signed_distance_mod(x_star, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
							vy_star = signed_distance_mod(y_star, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
							vz_star = signed_distance_mod(z_star, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
							(vx_star, vy_star, vz_star) = (	A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
															A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
															A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star)
							if abs(vx_star) <= R_ss[current_cell_list_ss[current_particle_in_cell], 1] && abs(vy_star) <= R_ss[current_cell_list_ss[current_particle_in_cell], 2] && abs(vz_star) <= R_ss[current_cell_list_ss[current_particle_in_cell], 3]
								is_proposed_position_ok = false
							end
						elseif particle_type_ss == "superellipsoid"
							vx_star = signed_distance_mod(x_star, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
							vy_star = signed_distance_mod(y_star, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
							vz_star = signed_distance_mod(z_star, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
							(vx_star, vy_star, vz_star) = (	A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
															A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
															A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star)
							if (abs(vx_star)/R_ss[current_cell_list_ss[current_particle_in_cell], 1])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] +
								(abs(vy_star)/R_ss[current_cell_list_ss[current_particle_in_cell], 2])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] +
								(abs(vz_star)/R_ss[current_cell_list_ss[current_particle_in_cell], 3])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] <= 1.0
								is_proposed_position_ok = false
							end
						end
					end

					# Make displacement if proposed position ok.
					if is_proposed_position_ok
						x = x_star
						y = y_star
						z = z_star

						x_abs = x_abs + deltax
						y_abs = y_abs + deltay
						z_abs = z_abs + deltaz

						D_empirical = D_empirical + deltax^2 + deltay^2 + deltaz^2

					end
				elseif boundary_condition == "multiple-rejection"
					is_proposed_position_ok = false
					while !is_proposed_position_ok
						deltax = sigma * randn()
						deltay = sigma * randn()
						deltaz = sigma * randn()

						x_star = position_mod(x + deltax, Lx)
						y_star = position_mod(y + deltay, Ly)
						z_star = position_mod(z + deltaz, Lz)

						# Check for diffuser-particle intersections.
						current_particle_in_cell = 0
						is_proposed_position_ok = true
						while current_particle_in_cell < number_of_particles_current_cell && is_proposed_position_ok
							current_particle_in_cell += 1

							if particle_type_ss == "sphere"
								vx_star = signed_distance_mod(x_star, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
								vy_star = signed_distance_mod(y_star, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
								vz_star = signed_distance_mod(z_star, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)

								if vx_star^2 + vy_star^2 + vz_star^2 <= R_ss[current_cell_list_ss[current_particle_in_cell], 1]^2
									is_proposed_position_ok = false
								end
							elseif particle_type_ss == "ellipse"
								# Make sure we can detect intersection by ensuring that v and v_star are close, not just close mod L.
							   	# Ellipse is a bit special in this regard.

								# Observe that the arguments to signed_distance_mod are shifted. Calling them
								# the same way as in other place will lead to wrong sign of distances.
								vx = signed_distance_mod(X_ss[current_cell_list_ss[current_particle_in_cell]], x, Lx)
								vy = signed_distance_mod(Y_ss[current_cell_list_ss[current_particle_in_cell]], y, Ly)
								vz = signed_distance_mod(Z_ss[current_cell_list_ss[current_particle_in_cell]], z, Lz)
								vx_star = vx + deltax
								vy_star = vy + deltay
								vz_star = vz + deltaz
								(vx, vy, vz) = (				A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx +	A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz,
													   			A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz,
													   			A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz)
								(vx_star, vy_star, vz_star) = (	A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
													   			A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
													   			A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star)

								alpha = - vz_star / (vz - vz_star)
								if 0.0 < alpha < 1.0
									w1_intersection = alpha * vx + (1.0 - alpha) * vx_star
									w2_intersection = alpha * vy + (1.0 - alpha) * vy_star
									if (w1_intersection/R_ss[current_cell_list_ss[current_particle_in_cell], 1])^2 + (w2_intersection/R_ss[current_cell_list_ss[current_particle_in_cell], 2])^2 < 1.0
										is_proposed_position_ok = false
									end
								end
							elseif particle_type_ss == "ellipsoid"
								vx_star = signed_distance_mod(x_star, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
								vy_star = signed_distance_mod(y_star, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
								vz_star = signed_distance_mod(z_star, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
								if vx_star * (A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star) + vy_star * (A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star) + vz_star * (A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star) <= 1.0
									is_proposed_position_ok = false
								end
							elseif particle_type_ss == "cuboid"
								vx_star = signed_distance_mod(x_star, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
								vy_star = signed_distance_mod(y_star, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
								vz_star = signed_distance_mod(z_star, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
								(vx_star, vy_star, vz_star) = (	A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
																A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
																A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star)
								if abs(vx_star) <= R_ss[current_cell_list_ss[current_particle_in_cell], 1] && abs(vy_star) <= R_ss[current_cell_list_ss[current_particle_in_cell], 2] && abs(vz_star) <= R_ss[current_cell_list_ss[current_particle_in_cell], 3]
									is_proposed_position_ok = false
								end
							elseif particle_type_ss == "superellipsoid"
								vx_star = signed_distance_mod(x_star, X_ss[current_cell_list_ss[current_particle_in_cell]], Lx)
								vy_star = signed_distance_mod(y_star, Y_ss[current_cell_list_ss[current_particle_in_cell]], Ly)
								vz_star = signed_distance_mod(z_star, Z_ss[current_cell_list_ss[current_particle_in_cell]], Lz)
								(vx_star, vy_star, vz_star) = (	A11_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A12_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A13_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
																A21_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A22_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A23_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star,
																A31_ss[current_cell_list_ss[current_particle_in_cell]] * vx_star + A32_ss[current_cell_list_ss[current_particle_in_cell]] * vy_star + A33_ss[current_cell_list_ss[current_particle_in_cell]] * vz_star)
								if (abs(vx_star)/R_ss[current_cell_list_ss[current_particle_in_cell], 1])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] +
									(abs(vy_star)/R_ss[current_cell_list_ss[current_particle_in_cell], 2])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] +
									(abs(vz_star)/R_ss[current_cell_list_ss[current_particle_in_cell], 3])^R_ss[current_cell_list_ss[current_particle_in_cell], 4] <= 1.0
									is_proposed_position_ok = false
								end

							end
						end
					end

					x = x_star
					y = y_star
					z = z_star

					x_abs = x_abs + deltax
					y_abs = y_abs + deltay
					z_abs = z_abs + deltaz

					D_empirical = D_empirical + deltax^2 + deltay^2 + deltaz^2
				end
			end

			trajectory_x[current_time_coarse] = x_abs
			trajectory_y[current_time_coarse] = y_abs
			trajectory_z[current_time_coarse] = z_abs

			ssd[current_time_coarse] += (trajectory_x[current_time_coarse] - trajectory_x[1])^2 + (trajectory_y[current_time_coarse] - trajectory_y[1])^2 + (trajectory_z[current_time_coarse] - trajectory_z[1])^2
			ssd_x[current_time_coarse] += (trajectory_x[current_time_coarse] - trajectory_x[1])^2
			ssd_y[current_time_coarse] += (trajectory_y[current_time_coarse] - trajectory_y[1])^2
			ssd_z[current_time_coarse] += (trajectory_z[current_time_coarse] - trajectory_z[1])^2
		end
	end
	output::Array{Float64, 1} = vcat(ssd, ssd_x, ssd_y, ssd_z, D_empirical)
	#output::Array{Float64, 1} = vcat(ssd, trajectory_x, trajectory_y, trajectory_z, D_empirical)
	return output
end
