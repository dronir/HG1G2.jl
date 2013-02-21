module HG1G2

using Splines
using PiecewiseFunctions
using HG1G2conversions

export default_basis

# Radian and degree conversion functions
rad{T<:Real}(x::Union(T,Array{T})) = 0.0174532925199433.*x
deg{T<:Real}(x::Union(T,Array{T})) = 57.2957795130823.*x

# Constants
const COEF_G1_SMALL = 0.7527
const CONST_G1_SMALL = 0.06164
const COEF_G2_SMALL = 0.9529
const CONST_G2_SMALL = 0.02162
const COEF_G1_LARGE = -0.9612
const CONST_G1_LARGE = 0.6270
const COEF_G2_LARGE = -0.6125
const CONST_G2_LARGE = 0.5572


function default_basis()
    basis_functions = Array(PiecewiseFunction, 3)
    
    # Construct linear part of a1
    a1_linear(x::Real) = -1.90985931710274*x + 1.0
    
    # Construct spline part of a1
    xvalues = rad([7.5, 30, 60, 90, 120, 150])
    yvalues = [0.75,0.33486,0.134106,0.0511048,0.0214657,0.0036397]
    deriv = [-1.90986, -0.0913286]
    S1 = Spline(xvalues, yvalues, deriv)
    a1_spline(x::Real) = SplineFunction(x, S1)

    basis_functions[1] = PiecewiseFunction()
    add_component!(basis_functions[1], a1_linear, 0.0, rad(7.5))
    add_component!(basis_functions[1], a1_spline, rad(7.5), rad(150))
    
    # Construct linear part of a2
    a2_linear(x::Real) = -0.572957795130823*x + 1.0
    
    # Construct spline part of a2
    xvalues = rad([7.5, 30, 60, 90, 120, 150])
    yvalues = [0.925,0.628842,0.317555,0.127164,0.0223739,0.000165057]
    deriv = [-0.572958, -8.6573138e-8]
    S2 = Spline(xvalues, yvalues, deriv)
    a2_spline(x::Real) = SplineFunction(x, S2)

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
    a3_spline(x::Real) = SplineFunction(x, S3)
    
    basis_functions[3] = PiecewiseFunction()
    add_component!(basis_functions[3], a3_constant, rad(30), rad(150))
    add_component!(basis_functions[3], a3_spline, 0.0, rad(30))
    return basis_functions
end




end #module
