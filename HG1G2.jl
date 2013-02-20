module HG1G2

using Splines
using PiecewiseFunctions
using HG1G2conversions

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







end #module
