function axis_aligned_bounding_box(	particle_type::String,
									R::Array{Float64, 1},
									x::Float64,
									y::Float64,
									z::Float64,
									q0::Float64,
									q1::Float64,
									q2::Float64,
									q3::Float64)
	lbx::Float64 = 0.0
	ubx::Float64 = 0.0
	lby::Float64 = 0.0
	uby::Float64 = 0.0
	lbz::Float64 = 0.0
	ubz::Float64 = 0.0

	if particle_type == "sphere"
		lbx = x - R[1]
		ubx = x + R[1]
		lby = y - R[1]
		uby = y + R[1]
		lbz = z - R[1]
		ubz = z + R[1]
	elseif particle_type == "ellipse"
		# Rotation matrix.
		(a11::Float64, a12::Float64, a13::Float64, a21::Float64, a22::Float64, a23::Float64, a31::Float64, a32::Float64, a33::Float64) = rotation_matrix(q0, q1, q2, q3)

		# Unrotated corner coordinates of the ellipse (effectively 2-D) AABB. We don't subtract the 1 from the z coordinates because we only use four corners now.
		C::Array{Float64, 2} = 2.0 * [[0.0, 0.0, 0.0] [0.0, 1.0, 0.0] [1.0, 0.0, 0.0] [1.0, 1.0, 0.0]]
		C[1, :] += - 1.0
		C[2, :] += - 1.0

		# Unrotated but correctly scaled corner coordinates of the ellipse AABB.
		for i = 1:2
			C[i, :] *= R[i]
		end

		# Rotated corner coordinates, yielding the ellipse OBB.
		for i = 1:8
			C[:, i] = [a11 a12 a13 ; a21 a22 a23 ; a31 a32 a33] * C[:, i]
		end

		# Lower and upper bounds of OBB.
		lbx = x + minimum(C[1, :])
		ubx = x + maximum(C[1, :])
		lby = y + minimum(C[2, :])
		uby = y + maximum(C[2, :])
		lbz = z + minimum(C[3, :])
		ubz = z + maximum(C[3, :])
	elseif particle_type == "ellipsoid"
		# Rotation matrix.
		(a11::Float64, a12::Float64, a13::Float64, a21::Float64, a22::Float64, a23::Float64, a31::Float64, a32::Float64, a33::Float64) = rotation_matrix(q0, q1, q2, q3)

		# Unrotated corner coordinates of the ellipsoids AABB.
		C::Array{Float64, 2} = 2.0 * [[0.0, 0.0, 0.0] [0.0, 0.0, 1.0] [0.0, 1.0, 0.0] [0.0, 1.0, 1.0] [1.0, 0.0, 0.0] [1.0, 0.0, 1.0] [1.0, 1.0, 0.0] [1.0, 1.0, 1.0]] - 1.0

		# Unrotated but correctly scaled corner coordinates of the ellipsoids AABB.
		for i = 1:3
			C[i, :] *= R[i]
		end

		# Rotated corner coordinates, yielding the ellipsoids OBB.
		for i = 1:8
			C[:, i] = [a11 a12 a13 ; a21 a22 a23 ; a31 a32 a33] * C[:, i]
		end

		# Lower and upper bounds of OBB.
		lbx = x + minimum(C[1, :])
		ubx = x + maximum(C[1, :])
		lby = y + minimum(C[2, :])
		uby = y + maximum(C[2, :])
		lbz = z + minimum(C[3, :])
		ubz = z + maximum(C[3, :])
	elseif particle_type == "cuboid"
		# Rotation matrix.
		(a11::Float64, a12::Float64, a13::Float64, a21::Float64, a22::Float64, a23::Float64, a31::Float64, a32::Float64, a33::Float64) = rotation_matrix(q0, q1, q2, q3)

		# Unrotated corner coordinates.
		C::Array{Float64, 2} = 2.0 * [[0.0, 0.0, 0.0] [0.0, 0.0, 1.0] [0.0, 1.0, 0.0] [0.0, 1.0, 1.0] [1.0, 0.0, 0.0] [1.0, 0.0, 1.0] [1.0, 1.0, 0.0] [1.0, 1.0, 1.0]] - 1.0

		# Unrotated but correctly scaled corner coordinates.
		for i = 1:3
			C[i, :] *= R[i]
		end

		# Rotated corner coordinates.
		for i = 1:8
			C[:, i] = [a11 a12 a13 ; a21 a22 a23 ; a31 a32 a33] * C[:, i]
		end

		# Lower and upper bounds of AABB.
		lbx = x + minimum(C[1, :])
		ubx = x + maximum(C[1, :])
		lby = y + minimum(C[2, :])
		uby = y + maximum(C[2, :])
		lbz = z + minimum(C[3, :])
		ubz = z + maximum(C[3, :])
	end

	return (lbx, ubx, lby, uby, lbz, ubz)
end
