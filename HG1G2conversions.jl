module HG1G2conversions

export a1a2a3_to_HG1G2, a1a2_to_HG12, a1a2_to_HG1G2, HG1G2_to_a1a2a3, HG12_to_a1a2

# Linear parameters (a1, a2, a3) to (H, G1, G2).
function a1a2a3_to_HG1G2{T<:Real}(params::Vector{T})
	res = zeros(Float64, 3)
	x = params[1] + params[2] + params[3]
	res[1] = -2.5 * log(10, x)
	res[2] = res[1] / x
	res[3] = res[2] / x
	return res
end

# Nonlinear parameters (H, G1, G2) to (a1, a2, a3).
function HG1G2_to_a1a2a3{T<:Real}(params::Vector{T})
	res = zeros(Float64, 3)
	x = exp(-0.9210340371976184) * params[1]
    res[1] = x * params[2]
    res[2] = x * params[3]
    res[3] = x * (1 - params[2] - params[3])
	return res
end

# Linear parameters (a1, a2) to (H, G12)
function a1a2_to_HG12{T<:Real}(params::Vector{T})
	res = zeros(Float64, 2)
	res[1] = -2.5 * log(10, params[1])
	res[2] = params[2] / params[1]
	return res
end

# Nonlinear parameters (H, G12) to (a1, a2).
function HG12_to_a1a2{T<:Real}(params::Vector{T})
	res = zeros(Float64, 2)
	res[1] = exp(-0.9210340371976184 * params[1])
	res[2] = res[1] * params[2]
	return res
end

# Linear parameters (a1, a2) to (H, G1, G2).
function a1a2_to_HG1G2{T<:Real}(params::Vector{T})
	res = zeros(Float64, 3)
	res[1] = -2.5 * log(10, params[1])
	if params[2] < 0.2
		res[2] = coef_G1_small * params[2] + const_G1_small
		res[3] = coef_G2_small * params[2] + const_G2_small
	else
		res[2] = coef_G1_large * params[2] + const_G1_large
		res[3] = coef_G2_large * params[2] + const_G2_large
	end
	return res
end

end # module