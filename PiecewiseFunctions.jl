module PiecewiseFunctions

export PiecewiseFunction
export add_component!, get_value

type PiecewiseFunction
    n_components::Int64
    limits::Matrix{Float64}
    functions::Array{Function}
end

# Empty constructor
PiecewiseFunction() = PiecewiseFunction(0, zeros(1,2), Function[])

# Push a new component to piecewise function
# TODO: check overlaps and continuity
function add_component!(F::PiecewiseFunction, f::Function, low::Real, hi::Real)
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
function get_value(F::PiecewiseFunction, x::Real)
    if F.n_components == 0
        error("Piecewise function has no components defined yet!")
    end
    i = 0
    for j = 1:F.n_components
        if F.limits[j,1] <= x <= F.limits[j,2]
            i = j
            break
        end
    end
    if i==0
        println("Warning: Value $x outside ranges of piecewise function!") 
        return NaN
    end
    return F.functions[i](x)
end

function get_value{T<:Real}(F::PiecewiseFunction, x::Array{T})
    n = size(x,1)
    result = zeros(n)
    for i = 1:n
        result[i] = get_value(F, x[i])
    end
    return result
end

end # module
