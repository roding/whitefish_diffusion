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

	elseif particle_type == "ellipsoid"

	elseif particle_type == "cuboid"

	end

	return (lbx, ubx, lby, uby, lbz, ubz)
end
