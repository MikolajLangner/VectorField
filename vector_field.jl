# module VectorField


using Makie
using LinearAlgebra
using DifferentialEquations
using ForwardDiff


# Supertype of ColorVectors (vectors with their colors, anchor points, etc).
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
function field2D(Fx::Function, Fy::Function;
                xbounds::Tuple{Real, Real} = (-1, 1),
                ybounds::Tuple{Real, Real} = (-1, 1))

    # Create points in ℝ²
    points = [Point2f0(i, j) for i in LinRange(xbounds[1], xbounds[2], 20)
                             for j in LinRange(ybounds[1], ybounds[2], 20)]

    # Define new function; `f` that takes points in ℝ²
    pointf(p::Point2f0) = [Fx(p...), Fy(p...)]

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
field2D((x,y)->y, (x, y)->x,
        xbounds = (-1,1), ybounds = (-1,1))


# Main fuction for 3D plotting
function field3D(Fx::Function, Fy::Function, Fz::Function;
                    xbounds::Tuple{Real, Real} = (-1, 1),
                    ybounds::Tuple{Real, Real} = (-1, 1),
                    zbounds::Tuple{Real, Real} = (-1, 1))

    # Create points in ℝ³
    points = [Point3f0(i, j, k) for i in LinRange(xbounds[1], xbounds[2], 10)
                                for j in LinRange(ybounds[1], ybounds[2], 10)
                                for k in LinRange(zbounds[1], zbounds[2], 10)]

    # Define new function; `f` that takes points in ℝ³
    pointf(p::Point3f0) = [Fx(p...), Fy(p...), Fz(p...)]

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
field3D((x,y,z)->y*z, (x,y,z)->x*z, (x, y, z)->x*y,
        xbounds = (-2,2), ybounds = (-2,2), zbounds = (-2,2))


# function for calculating position (x(t), y(t)) of an object
function positions2D(Fx::Function, Fy::Function;
                    startPoint::Array{T, 1} = [0, 0], # starting position
                    time::Tuple{Real, Real} = (0, 1), # time
                    timePoints::Integer = 100) where T <: Real # timepoints for interpolation

    # Define ∂x/∂t and ∂y/∂t
    function diffeqs(du, u, p, t)
        x, y = u
        du[1] = dx = Fx(x, y)
        du[2] = dy = Fy(x, y)
    end

    prob = ODEProblem(diffeqs, Float64.(startPoint), Float64.(time))
    sol = solve(prob)
    ts = LinRange(time[1], time[2], timePoints)
    X, Y = sol(ts, idxs = 1), sol(ts, idxs = 2)
    return X, Y # X - Array of x coords, Y - Array of y coords

end

# Example
X, Y = positions2D((x, y)->y, (x, y)->-x, startPoint=[-0.5, 0.5], time=(0, 5))
xk = X[end] # Current x position
yk = Y[end] # Current y position


# function for plotting trajectory based on calculations from "position2D"
function trajectory2D(Fx::Function, Fy::Function;
                        showField::Bool = :false,
                        startPoint::Array{T, 1} = [0, 0], # starting position
                        time::Tuple{Real, Real} = (0, 1), # time
                        timePoints::Integer = 100,
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        linewidth::Real = 3) where T <: Real

    X, Y = positions2D(Fx, Fy, startPoint = startPoint, time = time, timePoints = timePoints)
    X = X[xbounds[1] .< X .< xbounds[2]]
    Y = Y[ybounds[1] .< Y .< ybounds[2]]

    # X and Y must have the same length for plotting
    cut = min(length(X), length(Y))
    X, Y = X[1:cut], Y[1:cut]
    ts = LinRange(0, cut, 10*cut) # range only for color gradient

    if showField
        scene = field2D(Fx, Fy, xbounds = xbounds, ybounds = ybounds)
        lines!(scene, X, Y, linewidth=linewidth, color=ts, colormap=:grayC)
    else
        scene = lines(X, Y, linewidth=linewidth, color=ts, colormap=:grayC)
    end

    scatter!(scene, [X[end]], [Y[end]], markersize=xbounds[2]/20)

    scene |> display
    return scene

end

# Example

trajectory2D((x,y)->-y, (x,y)->x, xbounds=(-2,2), ybounds=(-2,2), linewidth=3, startPoint = [1, 1])


# function for calculating position (x(t), y(t), z(t)) of an object
function positions3D(Fx::Function, Fy::Function, Fz::Function;
                    startPoint::Array{T, 1} = [0, 0, 0], # starting position
                    time::Tuple{Real, Real} = (0, 1), # time
                    timePoints::Integer = 100) where T <: Real # timepoints for interpolation

    # Define ∂x/∂t, ∂y/∂t and ∂z/∂t
    function diffeqs(du, u, p, t)
        x, y, z = u
        du[1] = dx = Fx(x, y, z)
        du[2] = dy = Fy(x, y, z)
        du[3] = dz = Fz(x, y, z)
    end

    prob = ODEProblem(diffeqs, Float64.(startPoint), Float64.(time))
    sol = solve(prob)
    ts = LinRange(time[1], time[2], timePoints)
    X, Y, Z = sol(ts, idxs = 1), sol(ts, idxs = 2), sol(ts, idxs = 3)
    return X, Y, Z # X - Array of x coords, Y, Z similarly

end


# Example
X, Y, Z = positions3D((x,y,z)->2*x,(x,y,z)->-2*y,(x,y,z)->-2*z,
                        startPoint = [0.5, 0.5, -2], time=(0, 10), timePoints=500)
xk2 = X[end] # Current x position
yk2 = Y[end] # Current y position
zk2 = Z[end] # Current z position


# function for plotting trajectory based on calculations from "position3D"
function trajectory3D(Fx::Function, Fy::Function, Fz::Function;
                        showField::Bool = :false,
                        startPoint::Array{T, 1} = [0, 0, 0],
                        time::Tuple{Real, Real} = (0, 1),
                        timePoints::Integer = 100,
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        zbounds::Tuple{Real, Real} = (-1, 1),
                        linewidth::Real = 3) where T <: Real

    X, Y, Z = positions3D(Fx, Fy, Fz, startPoint = startPoint, time = time, timePoints = timePoints)

    X = X[xbounds[1] .< X .< xbounds[2]]
    Y = Y[ybounds[1] .< Y .< ybounds[2]]
    Z = Z[zbounds[1] .< Z .< zbounds[2]]

    # X, Y and Z must have the same length for plotting
    cut = min(length(X), length(Y), length(Z))
    X, Y, Z = X[1:cut], Y[1:cut], Z[1:cut]
    ts = LinRange(0, cut, 10*cut) # range only for color gradient

    if showField
        scene = field3D(Fx, Fy, Fz,
                    xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)
        lines!(scene, X, Y, Z, linewidth=linewidth, color=ts, colormap=:grayC)
    else
        scene = lines(X, Y, Z, linewidth=linewidth, color=ts, colormap=:grayC)
    end

    scatter!(scene, [X[end]], [Y[end]], [Z[end]], markersize=xbounds[2]/20)

    scene |> display
    return scene

end


function trajectory3D(X::RecursiveArrayTools.AbstractDiffEqArray,
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
    lines(X, Y, Z, linewidth=linewidth, color=ts, colormap=:grayC)
    meshscatter!([X[end]], [Y[end]], [Z[end]], markersize=xbounds[2]/30) |> display
    return AbstractPlotting.current_scene()
end
# Example
trajectory3D((x,y,z)->2*x,(x,y,z)->-2*y,(x,y,z)->-2*z,
                        startPoint = [0.5, 0.5, -2], time=(0, 100),
                        timePoints=1000, xbounds=(-2,2),ybounds=(-2,2),zbounds=(-2,2), linewidth=4)


# Plot a vector field from gradient of given two-variables function
function gradientField2D(f::Function,
                            xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real},
                            showContour::Bool=:true)

    # Create points from given intervals
    points = [Point2f0(i, j) for i in LinRange(xbounds[1], xbounds[2], 20)
                             for j in LinRange(ybounds[1], ybounds[2], 20)]

    # Define functions to calculate gradient
    vectorf(v::Vector) = f(v...)
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
# gradientField2D((x, y) -> x*exp(-x^2-y^2), (-2, 2), (-2, 2))


function gradientField3D(f::Function;
                            xbounds::Tuple{Real, Real} = (-1, 1),
                            ybounds::Tuple{Real, Real} = (-1, 1),
                            zbounds::Tuple{Real, Real} = (-1, 1),
                            showContour::Bool = :false)

    # Create points from given intervals
    points = [Point3f0(i, j, k) for i in LinRange(xbounds[1], xbounds[2], 10)
                                for j in LinRange(ybounds[1], ybounds[2], 10)
                                for k in LinRange(zbounds[1], zbounds[2], 10)]

    # Define functions to calculate gradient
    vectorf(v::Vector) = f(v...)
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

# end

# Example
# gradientField3D((x, y, z)->sin(x*y*z), (-1, 1), (-1, 1), (-1, 1))


function divergence(position::Array{T, 1}, fx::Function, fy::Function,
                    fz=missing) where T <: Real

    vectorFunctions = []
    vectorFuncX(v::Vector) = fx(v...)
    vectorFuncY(v::Vector) = fy(v...)
    push!(vectorFunctions, vectorFuncX)
    push!(vectorFunctions, vectorFuncY)
    if typeof(fz) <: Function
        vectorFuncZ(v::Vector) = fz(v...)
        push!(vectorFunctions, vectorFuncZ)
    end

    return sum([ForwardDiff.gradient(func, position)[variable]
                for (variable, func) in enumerate(vectorFunctions)])

end


function curl(position::Array{T, 1}, Fx::Function, Fy::Function,
                fz=missing) where T <: Real

    vectorFuncX(v::Vector) = Fx(v...)
    vectorFuncY(v::Vector) = Fy(v...)
    if typeof(Fz) <: Function
        vectorFuncZ(v::Vector) = Fz(v...)
        return [ForwardDiff.gradient(vectorFuncZ, position)[2] - ForwardDiff.gradient(vectorFuncY, position)[3],
                ForwardDiff.gradient(vectorFuncX, position)[3] - ForwardDiff.gradient(vectorFuncZ, position)[1],
                ForwardDiff.gradient(vectorFuncY, position)[1] - ForwardDiff.gradient(vectorFuncX, position)[2]]
    else
        return [0, 0, ForwardDiff.gradient(vectorFuncY, position)[1] - ForwardDiff.gradient(vectorFuncX, position)[2]]
    end

end


function animate2D(Fx::Function, Fy::Function, title::String;
                        showField::Bool = :false,
                        startPoints::Array = [[0, 0]], # starting position
                        time::Tuple{Real, Real} = (0, 1), # time
                        timePoints::Integer = 100,
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        fps::Integer = 24,
                        linewidth::Real = 3) where T <: Real

    xys = [] # Array for tuples: [(X1, Y1), (X2, Y2), ...]
    for startPoint in startPoints
        X, Y = positions2D(Fx, Fy, startPoint = startPoint, time = time, timePoints = timePoints)
        X = X[xbounds[1] .< X .< xbounds[2]]
        Y = Y[ybounds[1] .< Y .< ybounds[2]]
        cut = min(length(X), length(Y)) # X and Y must have the same length for plotting
        X, Y = X[1:cut], Y[1:cut] # X and Y must have the same length for plotting
        push!(xys, (X, Y)) # Add particle's coordinates to main Array
    end

    if showField
        scene = field2D(Fx, Fy, xbounds = xbounds, ybounds = ybounds)
        for xy in xys # some plots for future removing
            lines!([1], [1], visible=false)
            lines!([1], [1], visible=false)
        end
    else
        scene = lines(xbounds[1]:xbounds[2], ybounds[1]:ybounds[2])
        for xy in xys # some plots for future removing
            lines!(xbounds[1]:xbounds[2], ybounds[1]:ybounds[2])
        end
    end

    color_g = min(length(xys[1][1]), length(xys[1][2])) # color gradient
    ts = LinRange(0, color_g, 10*color_g) # range only for color gradient

    record(scene, title, 1:length(xys[1][1])-1; framerate = fps) do i
        for xy in xys # first removing: useless plots, every next time: old trajectory and dot
            delete!(scene, scene[end])
            delete!(scene, scene[end])
        end
        for xy in xys # draw new trajectory and new dot
          lines!(scene, xy[1][1:i], xy[2][1:i], linewidth=linewidth, color=ts, colormap=:grayC)
          scatter!(scene, [xy[1][i]], [xy[2][i]], markersize=xbounds[2]/20)
        end
    end
end

animate2DTEST((x, y)->y, (x, y)->-x, "Example2D.gif",
xbounds=(-1, 1), ybounds=(-1, 1), startPoints=[[0.5, 0.5], [-0.5, 0.5], [0.6, 0.2], [-0.4, -0.4], [0.1, -0.1]], time=(0., 5.), showField=:true)


function animate3D(
                   X::RecursiveArrayTools.AbstractDiffEqArray,
                   Y::RecursiveArrayTools.AbstractDiffEqArray,
                   Z::RecursiveArrayTools.AbstractDiffEqArray,
                   title::String;
                   xbounds::Tuple{Real, Real} = (-1, 1),
                   ybounds::Tuple{Real, Real} = (-1, 1),
                   zbounds::Tuple{Real, Real} = (-1, 1),
                   linewidth::Real = 3,
                   fps::Integer=24,
                   drawField::Bool=false,
                   Field::Scene=Scene())

    if drawField == true
      scene = Field
      lines!(xbounds[1]:xbounds[2], ybounds[1]:ybounds[2], zbounds[1]:zbounds[2])
      lines!(xbounds[1]:xbounds[2], ybounds[1]:ybounds[2], zbounds[1]:zbounds[2])
   else
       scene = lines(xbounds[1]:xbounds[2], ybounds[1]:ybounds[2], zbounds[1]:zbounds[2], visible=false)
       lines!(xbounds[1]:xbounds[2], ybounds[1]:ybounds[2], zbounds[1]:zbounds[2])
       lines!(xbounds[1]:xbounds[2], ybounds[1]:ybounds[2], zbounds[1]:zbounds[2])
   end
    record(scene, title, 1:length(X)-1; framerate = fps) do i
       delete!(scene, scene[end])
       delete!(scene, scene[end])
       trajectory3D!(X[1:i], Y[1:i], Z[1:i], xbounds=xbounds, ybounds=ybounds, zbounds=zbounds, linewidth=linewidth)

   end
end

X, Y, Z = position3D((x, y, z)-> y, (x, y, z)-> -x, (x, y, z)-> 2, xyz0=[0.5, -0.5, -1.5], tspan=(0., 7.))
animate3D(X, Y, Z, "Example3D.gif", xbounds=(-5, 5), ybounds=(-5, 5), zbounds=(-5, 5), fps=24)

# end
