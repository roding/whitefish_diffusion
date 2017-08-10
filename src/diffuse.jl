function diffuse(particle_type::String,
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
	D0::Float64,
	deltat_coarse::Float64,
	number_of_time_points_coarse::Int64,
	number_of_time_points_fine_per_coarse::Int64,
	number_of_diffusers::Int64)

	# Number of particles.
	number_of_particles::Int64 = length(X)

	# Standard deviation of Gaussian jumps.
	deltat_fine::Float64 = deltat_coarse / convert(Float64, number_of_time_points_fine_per_coarse)
	sigma::Float64 = sqrt(2 * D0 * deltat_fine)

	current_particle::Int64 = 0
	is_initial_position_ok::Bool = true
	is_proposed_position_ok::Bool = true

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
	trajectory_x::Array{Float64} = zeros(number_of_time_points_coarse)
	trajectory_y::Array{Float64} = zeros(number_of_time_points_coarse)
	trajectory_z::Array{Float64} = zeros(number_of_time_points_coarse)
	ssd::Array{Float64} = zeros(number_of_time_points_coarse)
	ssd_x::Array{Float64} = zeros(number_of_time_points_coarse)
	ssd_y::Array{Float64} = zeros(number_of_time_points_coarse)
	ssd_z::Array{Float64} = zeros(number_of_time_points_coarse)
	D0_empirical::Float64 = 0.0

	t_start_diffusion::Float64 = convert(Float64, time_ns()) / 1e9
	t_elapsed_diffusion::Float64 = 0.0

	chunk::Int64 = 0
	for current_diffuser = 1:number_of_diffusers
		if mod(current_diffuser, 100) == 0
			println(current_diffuser)
		end

		# Pick random, valid initial position.
		if particle_type == "sphere"
			is_initial_position_ok = false
			while !is_initial_position_ok
				x = Lx * rand()
				y = Ly * rand()
				z = Lz * rand()

				is_initial_position_ok = true
				current_particle = 0
				while current_particle < number_of_particles && is_initial_position_ok
					current_particle += 1
					vx = signed_distance_mod(x, X[current_particle], Lx)
					vy = signed_distance_mod(y, Y[current_particle], Lx)
					vz = signed_distance_mod(z, Z[current_particle], Lx)
					if vx^2 + vy^2 + vz^2 <= R[current_particle, 1]^2
						is_initial_position_ok = false
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
				# Random proposal displacement.
				deltax = sigma * randn()
				deltay = sigma * randn()
				deltaz = sigma * randn()

				# Calculate proposed new position.
				x_star = x + deltax
				y_star = y + deltay
				z_star = z + deltaz

				if x_star < 0.0
					x_star += Lx
				elseif x_star > Lx
					x_star -= Lx
				end
				if y_star < 0.0
					y_star += Ly
				elseif y_star > Ly
					y_star -= Ly
				end
				if z_star < 0.0
					z_star += Lz
				elseif z_star > Lz
					z_star -= Lz
				end

				# Check for diffuser-particle intersections.
				current_particle = 0
				is_proposed_position_ok = true
				while current_particle < number_of_particles && is_proposed_position_ok
					current_particle += 1

					if particle_type == "sphere"
						# Coordinates of candidate diffuser position relative to particle.
						vx_star = signed_distance_mod(x_star, X[current_particle], Lx)
						vy_star = signed_distance_mod(y_star, Y[current_particle], Ly)
						vz_star = signed_distance_mod(z_star, Z[current_particle], Lz)

						if vx_star^2 + vy_star^2 + vz_star^2 <= R[current_particle, 1]^2
							is_proposed_position_ok = false
						end
					end

					if is_proposed_position_ok
						x = x_star
						y = y_star
						z = z_star

						x_abs = x_abs + deltax
						y_abs = y_abs + deltay
						z_abs = z_abs + deltaz

						D0_empirical = D0_empirical + deltax^2 + deltay^2 + deltaz^2
					end
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

		#t_elapsed_diffusion = convert(Float64, time_ns()) / 1e9 - t_start_diffusion
		#if !silent_mode && convert(Int64, floor(t_elapsed_diffusion / 10.0)) > chunk
		#	print_progress(t_elapsed_diffusion, current_diffuser, number_of_diffusers)
		#end
		#chunk = convert(Int64, floor(t_elapsed_diffusion / 10.0))
	end

	output::Array{Float64, 1} = vcat(ssd, ssd_x, ssd_y, ssd_z, D0_empirical)
	return output
end
