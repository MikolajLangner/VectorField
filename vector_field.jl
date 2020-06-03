module VectorField


using Makie
using LinearAlgebra
using DifferentialEquations
using ForwardDiff


# Supertype of ColorVectors (vectors with their colors, anchorPoints, etc).
abstract type ColorVector end


# Type of 2D ColorVector
mutable struct ColorVector2D <: ColorVector

    magnitude::Float64  # Full magnitude, computed by `norm` function
    direction::Point2f0  # Representation of the vector, decreased to better visualization
    anchorPoint::Point2f0  # Point in ℝ², start of the vector
    color::RGBf0  # The color of the vector in RGB format (0 - black, 1 - white)
    size::Float64  # The proportionality of size to the longest vector's magnitude

    #= Inner constructor, vector is non decreased direction representation,
        decreaser tells how much it should be decreased, point is an anchorPoint =#
    ColorVector2D(vector::Array{T, 1}, point::Point2f0, decreaser::Real) where T <: Real = new(norm(vector),
                                                                        Point2f0(normalize(vector)/decreaser),
                                                                        point, RGBf0(0, 0, 0), 1.0)

end


mutable struct ColorVector3D <: ColorVector

    magnitude::Float64  # Full magnitude, computed by `norm` function
    direction::Point3f0  # Representation of the vector, decreased to better visualization
    anchorPoint::Point3f0  # Point in ℝ³, start of the vector
    color::RGBf0  # The color of the vector in RGB format (0 - black, 1 - white)
    size::Float64  # The proportionality of size to the longest vector's magnitude

    #= Inner constructor, vector is non decreased direction representation,
        decreaser tells how much it should be decreased, point is and anchorPoint =#
    ColorVector3D(vector::Array{T, 1}, point::Point3f0, decreaser::Real) where T <: Real = new(norm(vector),
                                                                        Point3f0(normalize(vector)/decreaser),
                                                                        point, RGBf0(0, 0, 0), 1.0)

end


# Overwrite base methods for ColorVector to make sort possible on it
Base.isless(v::ColorVector, w::ColorVector) = v.magnitude < w.magnitude ? true : false
Base.isequal(v::ColorVector, w::ColorVector) = v.magnitude == w.magnitude ? true : false


# Change the size of vector to the proportion of longest vector's magnitude
function changeLength(v::ColorVector, longest::Real)

    v.size = v.magnitude / longest
    v.direction = v.direction * v.size
    return v

end


#= Change the size of vector to the proportion of longest vector's magnitude,
    green - shortest vectors; red - longest, positive; blue - longest, negative =#
function changeColor(v::ColorVector, longest::Real)

    v.color = v.direction[end] >= 0 ? RGBf0(v.size, 1-v.size, 0) : RGBf0(0, 1-v.size, v.size)
    return v

end


# Change the vectors to indicate their values (proportional magnitudes)
function differVectors(vectors::Array{T, 1})::Array{T, 1} where T <: ColorVector

    vectors = sort(vectors)
    vectors = changeLength.(vectors, vectors[end].magnitude)
    vectors = changeColor.(vectors, vectors[end].magnitude)
    return vectors

end


# Add vector to a plot with its parameters
function plotVector(vector::ColorVector, scene::Scene)

    arrows!(scene, [vector.anchorPoint], [vector.direction],
            arrowsize=0.05vector.size, linewidth=4vector.size,
            linecolor=vector.color, arrowcolor=vector.color)
    return scene

end


# Main function for 2D plotting
function field2D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real})

    # Create points in ℝ²
    points = [Point2f0(i, j) for i in LinRange(xbounds[1], xbounds[2], 20)
                             for j in LinRange(ybounds[1], ybounds[2], 20)]

    # Define new function; `f` that takes points in ℝ²
    pointf(p::Point2f0) = f(p[1], p[2])

    # Create vectors
    vectors = @. ColorVector2D(pointf(points), points, 10)
    vectors = differVectors(vectors)

    # Create a plot
    scene = Scene()
    plotVector.(vectors, scene)

    # Display and return the plot
    scene |> display
    return scene

end


# Example
field2D((x,y)->[y,x],(-1,1),(-1,1))


# Main fuction for 3D plotting
function field3D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real}, zbounds::Tuple{Real,Real})

    # Create points in ℝ³
    points = [Point3f0(i, j, k) for i in LinRange(xbounds[1], xbounds[2], 10)
                                for j in LinRange(ybounds[1], ybounds[2], 10)
                                for k in LinRange(zbounds[1], zbounds[2], 10)]

    # Define new function; `f` that takes points in ℝ³
    pointf(p::Point3f0) = f(p[1], p[2], p[3])

    # Create vectors
    vectors = @. ColorVector3D(pointf(points), points, 2)
    vectors = differVectors(vectors)

    # Create a plot
    scene = Scene()
    plotVector.(vectors, scene)

    # Display and return the plot
    scene |> display
    return scene

end


# Example
field3D((x,y,z)->[y*z,x*z,x*y],(-2,2),(-2,2),(-2,2))


# function for calculating position (x(t), y(t)) of an object
function position2D(Vx::Function, Vy::Function;
                    xy0::Array{Float64, 1} = [0, 0], # starting position
                    tspan::Tuple{Float64, Float64} = (0, 1), # time
                    timepoints::Integer = 100) # timepoints for interpolation

    # Define ∂x/∂t and ∂y/∂t
    function diffeqs(du,u, p, t)
        x,y = u
        du[1] = dx = Vx(x, y)
        du[2] = dy = Vy(x, y)
    end

    prob = ODEProblem(diffeqs, xy0, tspan)
    sol = solve(prob)
    ts = LinRange(tspan[1],tspan[2], timepoints)
    X, Y = sol(ts,idxs=1), sol(ts,idxs=2)
    return X, Y # X - Array of x coords, Y - Array of y coords
end

# Example
X, Y = position2D((x, y)->y, (x, y)->-x, xy0=[-0.5, 0.5], tspan=(0., 5.))
xk = X[end] # Current x position
yk = Y[end] # Current y position


# function for plotting trajectory based on calculations from "position2D"
function trajectory2D!(X::RecursiveArrayTools.AbstractDiffEqArray,
                       Y::RecursiveArrayTools.AbstractDiffEqArray;
                       xbounds::Tuple{Real, Real} = (-1, 1),
                       ybounds::Tuple{Real, Real} = (-1, 1),
                       linewidth::Real = 3)

    X = X[xbounds[1] .< X .< xbounds[2]]
    Y = Y[ybounds[1] .< Y .< ybounds[2]]

    # X and Y must have the same length for plotting
    cut = min(length(X), length(Y))
    X, Y = X[1:cut], Y[1:cut]
    ts = LinRange(0, cut, 10*cut) # range only for color gradient
    scatter!([X[end]], [Y[end]], markersize=xbounds[2]/20)
    lines!(X, Y, linewidth=linewidth, color=ts, colormap=:grayC) |> display
end

# Example

field2D((x,y)->[y,-x],(-2,2),(-2,2))
trajectory2D!(X, Y, xbounds=(-2,2), ybounds=(-2,2), linewidth=3)


# function for calculating position (x(t), y(t), z(t)) of an object
function position3D(Vx::Function, Vy::Function, Vz::Function;
                    xyz0::Array{Float64, 1} = [0, 0, 0], # starting position
                    tspan::Tuple{Float64, Float64} = (0, 1), # time
                    timepoints::Integer = 100) # timepoints for interpolation

    # Define ∂x/∂t, ∂y/∂t and ∂z/∂t
    function diffeqs(du,u, p, t)
        x,y,z = u
        du[1] = dx = Vx(x, y, z)
        du[2] = dy = Vy(x, y, z)
        du[3] = dz = Vz(x, y, z)
    end

    prob = ODEProblem(diffeqs, xyz0,tspan)
    sol = solve(prob)
    ts = LinRange(tspan[1], tspan[2], timepoints)
    X, Y, Z = sol(ts,idxs=1), sol(ts,idxs=2), sol(ts,idxs=3)
    return X, Y, Z # X - Array of x coords, Y, Z similarly
end

# Example

X, Y, Z = position3D((x,y,z)->2*x,(x,y,z)->-2*y,(x,y,z)->-2*z,
                    xyz0 = [0.5, 0.5, -2.], tspan=(0., 10.0),
                    timepoints=500)
xk2 = X[end] # Current x position
yk2 = Y[end] # Current y position
zk2 = Z[end] # Current z position


# function for plotting trajectory based on calculations from "position3D"
function trajectory3D!(X::RecursiveArrayTools.AbstractDiffEqArray,
                       Y::RecursiveArrayTools.AbstractDiffEqArray,
                       Z::RecursiveArrayTools.AbstractDiffEqArray;
                       xbounds::Tuple{Real, Real} = (-1, 1),
                       ybounds::Tuple{Real, Real} = (-1, 1),
                       zbounds::Tuple{Real, Real} = (-1, 1),
                       linewidth::Real = 3)

    X = X[xbounds[1] .< X .< xbounds[2]]
    Y = Y[ybounds[1] .< Y .< ybounds[2]]
    Z = Z[zbounds[1] .< Z .< zbounds[2]]

    # X, Y and Z must have the same length for plotting
    cut = min(length(X), length(Y), length(Z))
    X, Y, Z = X[1:cut], Y[1:cut], Z[1:cut]
    ts = LinRange(0, cut, 10*cut) # range only for color gradient
    meshscatter!([X[end]], [Y[end]], [Z[end]], markersize=xbounds[2]/20)
    lines!(X, Y, Z, linewidth=linewidth, color=ts, colormap=:grayC) |> display
end

# Example

field3D((x,y,z)->[2*x,-2*y,-2*z], (-2,2),(-2,2),(-2,2))
trajectory3D!(X, Y, Z, xbounds=(-2,2),ybounds=(-2,2),zbounds=(-2,2), linewidth=4)


# Plot a vector field from gradient of given two-variables function
function gradientField2D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real},
                         showContour::Bool=:true)

    # Create points from given intervals
    points = [Point2f0(i, j) for i in LinRange(xbounds[1], xbounds[2], 20)
                             for j in LinRange(ybounds[1], ybounds[2], 20)]

    # Define functions to calculate gradient
    vectorf(v::Vector) = f(v[1], v[2])
    ∇(p::Point2f0) = ForwardDiff.gradient(vectorf, [p...])

    # Create vectors with specified anchor points
    vectors = @. ColorVector2D(∇(points), points, 10)
    vectors = differVectors(vectors)

    # Create a plot
    scene = Scene()
    plotVector.(vectors, scene)

    # Plot contour, if needed
    if showContour
        contour!(LinRange(xbounds[1], xbounds[2], 100),
                 LinRange(ybounds[1], ybounds[2], 100),
                 f.(LinRange(xbounds[1], xbounds[2], 100), LinRange(ybounds[1], ybounds[2], 100)'))
    end

    # Display and return the plot
    scene |> display
    return scene

end


# Example
gradientField2D((x, y) -> x*exp(-x^2-y^2), (-2, 2), (-2, 2))


function gradientField3D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real},
                         zbounds::Tuple{Real, Real}, showContour::Bool=:true)

    # Create points from given intervals
    points = [Point3f0(i, j, k) for i in LinRange(xbounds[1], xbounds[2], 10)
                                for j in LinRange(ybounds[1], ybounds[2], 10)
                                for k in LinRange(zbounds[1], zbounds[2], 10)]

    # Define functions to calculate gradient
    vectorf(v::Vector) = f(v[1], v[2], v[3])
    ∇(p::Point3f0) = ForwardDiff.gradient(vectorf, [p...])

    # Create vectors with specified anchor points
    vectors = @. ColorVector3D(∇(points), points, 10)
    vectors = differVectors(vectors)

    # Create a plot
    scene = Scene()
    plotVector.(vectors, scene)

    # Display and return the plot
    scene |> display
    return scene

end

# Example
gradientField3D((x, y, z)->sin(x*y*z), (-1, 1), (-1, 1), (-1, 1))

end
