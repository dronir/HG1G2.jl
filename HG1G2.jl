module HG1G2

export a1a2a3_to_HG1G2, a1a2_to_HG12, a1a2_to_HG1G2, HG1G2_to_a1a2a3, HG12_to_a1a2
export Spline, spline_coefs
export PieceWiseFunction, add_component!, get_value
export ConstantFunction, LinearFunction, SplineFunction

# Radian and degree conversion functions
rad{T<:Real}(x::Union(T,Array{T})) = 0.0174532925199433.*x
deg{T<:Real}(x::Union(T,Array{T})) = 57.2957795130823.*x

# Constants
const coef_G1_small = 0.7527
const const_G1_small = 0.06164
const coef_G2_small = 0.9529
const const_G2_small = 0.02162
const coef_G1_large = -0.9612
const const_G1_large = 0.6270
const coef_G2_large = -0.6125
const const_G2_large = 0.5572



ConstantFunction(x::Real, value::Real) = value
ConstantFunction{T<:Real}(x::Array{T}, value::Real) = value*ones(typeof(value),size(L))
LinearFunction{T<:Real}(x::Union(T,Array{T}), a::Real, b::Real) = a.*x + b

type Spline
	knots::Vector{Float64}
	coef::Matrix{Float64}
end

function SplineFunction(x::Real, spline::Spline)
	n_knots = size(spline.knots,1)
	j = 0
	for j = 1:n_knots-1
		if spline.knots[j] <= x <= spline.knots[j+1]
			break
		end
	end
	t = (x - spline.knots[j]) / (spline.knots[j+1] - spline.knots[j])
	return spline.coef[j,1] + spline.coef[j,2]*t + spline.coef[j,3]*t^2 + spline.coef[j,4]*t^3
end

type PieceWiseFunction
	n_components::Int64
	limits::Matrix{Float64}
	functions::Array{Function}
end

# Dummy constructor
PieceWiseFunction() = PieceWiseFunction(0, zeros(1,2), Function[])

# Push a new component to piecewise function
# TODO: check overlaps and continuity
function add_component!(F::PieceWiseFunction, f::Function, low::Real, hi::Real)
	if F.n_components == 0
		F.limits[1,1] = low
		F.limits[1,2] = hi
	else
		F.limits = vcat(F.limits, [low, hi]') #' syntax colour bug
	end
	push!(F.functions, f)
	F.n_components += 1
end

# Get value of piecewise function at x
function get_value(F::PieceWiseFunction, x::Real)
	if F.n_components == 0
		error("Piecewise function has no components defined yet!")
	end
	j = 0
	for j = 1:F.n_components
		if F.limits[j,1] <= x <= F.limits[j,2] break end
	end
	if j==0
		println("Warning: Value $x outside ranges of piecewise function!") 
		return NaN
	end
	return F.functions[j](x)
end

function get_value{T<:Real}(F::PieceWiseFunction, x::Array{T})
	n = size(x,1)
	result = zeros(n)
	for i = 1:n
		result[i] = get_value(F, x[i])
	end
	return result
end



function spline_coefs(xval::Vector{Float64}, yval::Vector{Float64}, deriv::Vector{Float64})
	N = size(xval, 1)
	
	A = ones(Float64, N)
	B = 4*ones(Float64, N)
	C = ones(Float64, N)
	A[1] = 0.0
	A[N] = 0.0
	B[1] = 1.0
	B[N] = 1.0
	C[1] = 0.0
	C[N] = 0.0

	R = zeros(Float64, N)
	gamma = zeros(Float64, N)
	U = zeros(Float64, N)	

	R[1] = deriv[1] * (xval[2] - xval[1])
	for i = 2:N-1
		R[i] = 3 * (yval[i+1] - yval[i-1])
	end
	R[N] = deriv[2] * (xval[N] - xval[N-1])
	
	# Trilinear solver
	beta = B[1]
	U[1] = R[1]/beta
	for i = 2:N
		gamma[i] = C[i-1] / beta
		beta = B[i] - A[i] * gamma[i]
		U[i] = (R[i] - A[i]*U[i-1]) / beta
	end
	for i = N-1 : -1 : 1
		U[i] = U[i] - gamma[i+1]*U[i+1]
	end
	
	S = Spline(xval, zeros(Float64, N-1, 4))
	for i = 1:N-1
		S.coef[i,1] = yval[i]
		S.coef[i,2] = U[i]
		S.coef[i,3] = 3 * (yval[i+1] - yval[i]) - 2 * U[i] - U[i+1]
		S.coef[i,4] = 2 * (yval[i] - yval[i+1]) + U[i] + U[i+1]
	end
	return S
end


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


end #module
