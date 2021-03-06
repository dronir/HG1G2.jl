module HG1G2

using Splines
using PiecewiseFunctions
using HG1G2conversions

export BasisFunctions
export form_base, fit_HG1G2, fitted_curve

# Radian and degree conversion functions
rad{T<:Real}(x::Union(T,Array{T})) = 0.017453292519943295.*x
deg{T<:Real}(x::Union(T,Array{T})) = 57.29577951308232.*x

type BasisFunctions
    funcs::Array{PiecewiseFunction}
    version::Int64
    lower::Float64
    upper::Float64
end



form_base(filename::String) = error("Not implemented yet.")

function form_base()
    basis_functions = Array(PiecewiseFunction, 3)
    
    # Construct linear part of a1
    a1_linear(x::Real) = 1.0 - 1.90985931710274*x
    
    # Construct spline part of a1
    xvalues = rad([7.5, 30, 60, 90, 120, 150])
    yvalues = [0.75, 0.33486, 0.134106, 0.0511048, 0.0214657, 0.0036397]
    deriv = [-1.90986, -0.0913286]
    S1 = Spline(xvalues, yvalues, deriv)
    a1_spline(x::Real) = SplineFunction(S1, x)

    basis_functions[1] = PiecewiseFunction()
    add_component!(basis_functions[1], a1_linear, 0.0, rad(7.5))
    add_component!(basis_functions[1], a1_spline, rad(7.5), rad(150))
    
    # Construct linear part of a2
    a2_linear(x::Real) = 1.0 - 0.572957795130823*x
    
    # Construct spline part of a2
    xvalues = rad([7.5, 30, 60, 90, 120, 150])
    yvalues = [0.925,0.628842,0.317555,0.127164,0.0223739,0.000165057]
    deriv = [-0.572958, -8.6573138e-8]
    S2 = Spline(xvalues, yvalues, deriv)
    a2_spline(x::Real) = SplineFunction(S2, x)

    basis_functions[2] = PiecewiseFunction()
    add_component!(basis_functions[2], a2_linear, 0.0, rad(7.5))
    add_component!(basis_functions[2], a2_spline, rad(7.5), rad(150))
    
    
    # Construct constant part of a3
    a3_constant(x::Real) = 0.0
    
    # Construct spline part of a3
    xvalues = rad([0.0, 0.3, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0, 30.0])
    yvalues = [1.,0.833812,0.577354,0.421448,0.231742,0.103482,0.0617335,0.016107,0.0]
    deriv = [-0.106301, 0.0]
    S3 = Spline(xvalues, yvalues, deriv)
    a3_spline(x::Real) = SplineFunction(S3, x)
        
    basis_functions[3] = PiecewiseFunction()
    add_component!(basis_functions[3], a3_spline, 0.0, rad(30))
    add_component!(basis_functions[3], a3_constant, rad(30), rad(150))
    BasisFunctions(basis_functions, 20101000, 0.0, rad(150))
end


function fit_HG1G2{T<:Real}(basis::BasisFunctions, data::Matrix{T}, errors::Vector{T})
    Ndata = size(data,1)
    xvalues = vec(data[:,1]) * pi/180
    yvalues = 10 .^(-0.4 * vec(data[:,2]))
    
    sigmas = yvalues .* (10.^(0.4 * errors) - 1)
    
    Nfuncs = size(basis.funcs,1)
    Amatrix = zeros(Ndata, Nfuncs)
    for i = 1:Ndata
        for j = 1:Nfuncs
            Amatrix[i,j] = get_value(basis.funcs[j], xvalues[i]) / sigmas[i]
        end
    end
    
    yvalues = 1/(10.^(0.4 * errors) - 1)

    as = Amatrix \ yvalues
    return a1a2a3_to_HG1G2(vec(as))
end

fit_HG1G2{T<:Real}(basis::BasisFunctions, data::Matrix{T}) = fit_HG1G2(basis,data,0.03*ones(size(data,1)))


function fitted_curve(T, params, basis)
    H = params[1]
    G1 = params[2]
    G2 = params[3]
    V0 = 10^(-0.4*H)
    Y1 = G1 * get_value(basis[1], T)
    Y2 = G2 * get_value(basis[2], T)
    Y3 = (1-G1-G2) * get_value(basis[3], T)
    return H - log(10, Y1 + Y2 + Y3) / 0.4
end


end #module
