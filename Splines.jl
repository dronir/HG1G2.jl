module Splines

export Spline, SplineFunction

type Spline
    knots::Vector{Float64}
    coef::Matrix{Float64}
end

function Spline(xval::Vector{Float64}, yval::Vector{Float64}, deriv::Vector{Float64})
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
    
    # tri-diagonal solver
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

end # module