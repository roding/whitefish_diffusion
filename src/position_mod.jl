function position_mod(x::Float64, L::Float64)
	if x < 0.0
		x += L
	elseif x > L
		x -= L
	end

	return x
end
