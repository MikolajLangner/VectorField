module VectorField


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
                                                                        Point2f0(normalize(vector) / decreaser),
                                                                        point, RGBf0(0, 0, 0), 1.0)

end


# Type of 3D ColorVector
mutable struct ColorVector3D <: ColorVector

    magnitude::Float64  # Full magnitude, computed by `norm` function
    direction::Point3f0  # Representation of the vector, decreased to better visualization
    anchorPoint::Point3f0  # Point in ℝ³, start of the vector
    color::RGBf0  # The color of the vector in RGB format (0 - black, 1 - white)
    size::Float64  # The proportionality of size to the longest vector's magnitude

    #= Inner constructor, vector is non decreased direction representation,
        decreaser tells how much it should be decreased, point is and anchorPoint =#
    ColorVector3D(vector::Array{T, 1}, point::Point3f0, decreaser::Real) where T <: Real = new(norm(vector),
                                                                        Point3f0(normalize(vector) / decreaser),
                                                                        point, RGBf0(0, 0, 0), 1.0)

end


# Overwrite base methods for ColorVector to make sort possible on it
Base.isless(v::ColorVector, w::ColorVector) = v.magnitude < w.magnitude ? true : false
Base.isequal(v::ColorVector, w::ColorVector) = v.magnitude == w.magnitude ? true : false


# Supertype for objects moving in a field
abstract type Body end


# Type of 2D object
mutable struct Body2D <: Body

    X::Array  # x coordinate positions
    Y::Array  # y coordinate positions
    colors::LinRange  # size of colors' gradient
    inBounds::Bool  # flag that indicates if object is in bounds at any time

    # Inner constructor for given tuple of arrays that represents positions
    Body2D(t::Tuple) = new(t[1], t[2], LinRange(0, length(t[1]), 10length(t[1])), :true)

end

# Type of 3D object
mutable struct Body3D <: Body

    X::Array  # x coordinate positions
    Y::Array  # y coordinate positions
    Z::Array  # z coordinate positions
    colors::LinRange  # size of colors' gradient
    inBounds::Bool  # flag that indicates if object is in bounds at any time

    # Inner constructor for given tuple of arrays that represents positions
    Body3D(t::Tuple) = new(t[1], t[2], t[3], LinRange(0, length(t[1]), 10length(t[1])), :true)

end


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
function plotVector(vector::ColorVector, scene::Scene;
                    xbounds::Tuple{Real, Real} = (-1, 1),
                    ybounds::Tuple{Real, Real} = (-1, 1),
                    zbounds::Tuple{Real, Real} = (0, 0))

    sumBounds = sum([bounds[2]-bounds[1] for bounds in [xbounds, ybounds, zbounds]])  # parameter for size of arrows
    arrows!(scene, [vector.anchorPoint], [vector.direction],
            arrowsize = 0.01vector.size*sumBounds, linewidth = 2vector.size,
            linecolor = vector.color, arrowcolor = vector.color)
    return scene

end


"""
    field(Fx::Function, Fy::Function; kwargs...)

Plot a vector field of two-variable functions.

# Arguments
- `xbounds::Tuple{Real, Real} = (-1, 1)`: boundaries for x coordinates.
- `ybounds::Tuple{Real, Real} = (-1, 1)`: boundaries for y coordinates.
"""
function field(Fx::Function, Fy::Function;
                xbounds::Tuple{Real, Real} = (-1, 1),
                ybounds::Tuple{Real, Real} = (-1, 1))

    # Create points in ℝ²
    points = [Point2f0(i, j) for i in LinRange(xbounds[1], xbounds[2], 20)
                             for j in LinRange(ybounds[1], ybounds[2], 20)]

    # Define new function; `f` that takes points in ℝ²
    pointf(p::Point2f0) = [Fx(p...), Fy(p...)]

    # Create vectors
    vectors = @. ColorVector2D(pointf(points), points, 40/((xbounds[2]-xbounds[1])+(ybounds[2]-ybounds[1])))
    vectors = differVectors(vectors)

    # Create a plot
    scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], visible = :false)
    plotVector.(vectors, scene, xbounds = xbounds, ybounds = ybounds)

    # Display and return the plot
    scene |> display
    return scene

end


"""
    field(Fx::Function, Fy::Function, Fz::Function; kwargs...)

Also for three-variable functions.

# Arguments
- `zbounds::Tuple{Real, Real} = (-1, 1)`: boundaries for z coordinates.
"""
function field(Fx::Function, Fy::Function, Fz::Function;
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
    sumBounds = sum([bounds[2]-bounds[1] for bounds in [xbounds, ybounds, zbounds]])
    vectors = @. ColorVector3D(pointf(points), points, 30/sumBounds)
    vectors = differVectors(vectors)

    # Create a plot
    scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], [zbounds[1], zbounds[2]], visible = :false)
    plotVector.(vectors, scene, xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)

    # Display and return the plot
    scene |> display
    return scene

end


# Keep only positions in boundaries
function cutPositions(body::Body;
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        zbounds::Tuple{Real, Real} = (-1, 1))

    # Check which points are in boundaries
    body.X = body.X[xbounds[1] .< body.X .< xbounds[2]]
    body.Y = body.Y[ybounds[1] .< body.Y .< ybounds[2]]
    if typeof(body) == Body3D
        body.Z = body.Z[zbounds[1] .< body.Z .< zbounds[2]]

        # X, Y and Z must have the same length for plotting
        cut = min(length(body.X), length(body.Y), length(body.Z))
        if cut > 0
            body.X = body.X[1:cut]
            body.Y = body.Y[1:cut]
            body.Z = body.Z[1:cut]
            body.colors = LinRange(0, cut, 10cut)
        else
            body.inBounds = :false
        end
    else
        cut = min(length(body.X), length(body.Y))
        if cut > 0
            body.X = body.X[1:cut]
            body.Y =  body.Y[1:cut]
            body.colors = LinRange(0, cut, 10cut)
        else
            body.inBounds = :false
        end
    end

    return body

end


# Plot objects with given specifications
function plotBody(body::Body, scene::Scene, linewidth::Real; stopFrame::Integer = length(body.X))

    # Check whether any of the points is in boundaries
    if body.inBounds
        if typeof(body) == Body3D
            lines!(scene, body.X[1:stopFrame], body.Y[1:stopFrame], body.Z[1:stopFrame],
                    linewidth = linewidth, color = body.colors[1:10stopFrame], colormap = :grayC)
            scatter!(scene, [body.X[stopFrame]], [body.Y[stopFrame]], [body.Z[stopFrame]], markersize=linewidth/70)
        else
            lines!(scene, body.X[1:stopFrame], body.Y[1:stopFrame],
                    linewidth = linewidth, color = body.colors[1:10stopFrame], colormap = :grayC)
            scatter!(scene, [body.X[stopFrame]], [body.Y[stopFrame]], markersize=linewidth/70)
        end
    end

    return scene

end


"""
    positions(Fx::Function, Fy::Function, startPoint::Array{T, 1} = [0, 0]; kwargs...) where T <: Real

Compute positions of an object for vector field of given two-variable functions.

# Arguments
- `time::Tuple{Real, Real} = (0, 1)`: boundaries for time which positions will be computed in.
- `timePoints::Integer = 100`: number of positions (accuracy) to compute.
"""
function positions(Fx::Function, Fy::Function, startPoint::Array{T, 1} = [0, 0]; # starting position
                    time::Tuple{Real, Real} = (0, 1), # time
                    timePoints::Integer = 100) where T <: Real # timepoints for interpolation

    # Define ∂x/∂t and ∂y/∂t
    function diffeqs(du, u, p, t)
        x, y = u
        du[1] = dx = Fx(x, y)
        du[2] = dy = Fy(x, y)
    end

    # Define the `problem` and solve it
    prob = ODEProblem(diffeqs, Float64.(startPoint), Float64.(time))
    sol = solve(prob)
    ts = LinRange(time[1], time[2], timePoints)
    X, Y = sol(ts, idxs = 1), sol(ts, idxs = 2)
    return X, Y # X - Array of x coords, Y - Array of y coords

end


"""
    positions(Fx::Function, Fy::Function, Fz::Function, startPoint::Array{T, 1} = [0, 0, 0]; kwargs...) where T <: Real

Also for three-variable functions.
"""
function positions(Fx::Function, Fy::Function, Fz::Function,
                    startPoint::Array{T, 1} = [0, 0, 0]; # starting position
                    time::Tuple{Real, Real} = (0, 1), # time
                    timePoints::Integer = 100) where T <: Real # timepoints for interpolation

    # Define ∂x/∂t, ∂y/∂t and ∂z/∂t
    function diffeqs(du, u, p, t)
        x, y, z = u
        du[1] = dx = Fx(x, y, z)
        du[2] = dy = Fy(x, y, z)
        du[3] = dz = Fz(x, y, z)
    end

    # Define the `problem` and solve it
    prob = ODEProblem(diffeqs, Float64.(startPoint), Float64.(time))
    sol = solve(prob)
    ts = LinRange(time[1], time[2], timePoints)
    X, Y, Z = sol(ts, idxs = 1), sol(ts, idxs = 2), sol(ts, idxs = 3)
    return X, Y, Z # X - Array of x coords, Y, Z similarly

end


"""
    trajectory(Fx::Function, Fy::Function, startPoints::Array{Array{T, 1}, 1}; kwargs...) where T <: Real

Plot full trajectories for objects with given `startPoints` and two-variable functions.

# Arguments
- `showField::Bool = :false`: flag to indicate whether plot vector field (arrows).
- `time::Tuple{Real, Real} = (0, 1)`: boundaries for time which positions will be computed in.
- `timePoints::Integer = 100`: number of positions (accuracy) to compute.
- `xbounds::Tuple{Real, Real} = (-1, 1)`: boundaries for x coordinates.
- `ybounds::Tuple{Real, Real} = (-1, 1)`: boundaries for y coordinates.
- `linewidth::Real = 3`: width of trajectories' lines.
"""
function trajectory(Fx::Function, Fy::Function,
                        startPoints::Array{Array{T, 1}, 1};
                        showField::Bool = :false,
                        time::Tuple{Real, Real} = (0, 1),
                        timePoints::Integer = 100,
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        linewidth::Real = 3) where T <: Real

    # Create objects and keep only positions in boundaries
    bodies = @. Body2D(positions(Fx, Fy, startPoints, time = time, timePoints = timePoints))
    bodies = cutPositions.(bodies, xbounds = xbounds, ybounds = ybounds)

    # Create a plot (with or without vector field)
    if showField
        scene = field(Fx, Fy, xbounds = xbounds, ybounds = ybounds)
    else
        scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], visible = :false)
    end

    # Plot objects
    plotBody.(bodies, scene, linewidth)

    # Display and return the plot
    scene |> display
    return scene

end


"""
    trajectory(Fx::Function, Fy::Function, Fz::Function, startPoints::Array{Array{T, 1}, 1}; kwargs...) where T <: Real

Also for three-variable functions.

# Arguments
- `zbounds::Tuple{Real, Real} = (-1, 1)`: boundaries for z coordinates.
"""
function trajectory(Fx::Function, Fy::Function, Fz::Function,
                        startPoints::Array{Array{T, 1}, 1};
                        showField::Bool = :false,
                        time::Tuple{Real, Real} = (0, 1),
                        timePoints::Integer = 100,
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        zbounds::Tuple{Real, Real} = (-1, 1),
                        linewidth::Real = 3) where T <: Real

    # Create objects and keep only positions in boundaries
    bodies = @. Body3D(positions(Fx, Fy, Fz, startPoints, time = time, timePoints = timePoints))
    bodies = cutPositions.(bodies, xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)

    # Create a plot (with or without vector field)
    if showField
        scene = field(Fx, Fy, Fz,
                    xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)
    else
        scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], [zbounds[1], zbounds[2]], visible = :false)
    end

    # Plot objects
    plotBody.(bodies, scene, linewidth)

    # Display and return the plot
    scene |> display
    return scene

end


function gradientField2D(f::Function; showContour::Bool = :false,
                            xbounds::Tuple{Real, Real} = (-1, 1),
                            ybounds::Tuple{Real, Real} = (-1, 1))

    # Create points from given intervals
    points = [Point2f0(i, j) for i in LinRange(xbounds[1], xbounds[2], 20)
                             for j in LinRange(ybounds[1], ybounds[2], 20)]

    # Define functions to calculate gradient
    vectorf(v::Vector) = f(v...)
    ∇(p::Point2f0) = ForwardDiff.gradient(vectorf, [p...])

    # Create vectors with specified anchor points
    vectors = @. ColorVector2D(∇(points), points, 40/((xbounds[2]-xbounds[1])+(ybounds[2]-ybounds[1])))
    vectors = differVectors(vectors)

    # Create a plot
    scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], visible = :false)
    plotVector.(vectors, scene, xbounds = xbounds, ybounds = ybounds)

    # Plot contour, if needed
    if showContour
        xs = LinRange(xbounds[1], xbounds[2], 100Integer(xbounds[2]-xbounds[1]))
        ys = LinRange(ybounds[1], ybounds[2], 100Integer(ybounds[2]-ybounds[1]))
        contour!(xs, ys, f.(xs, ys'))
    end

    # Display and return the plot
    scene |> display
    return scene

end


function gradientField3D(f::Function;
                            xbounds::Tuple{Real, Real} = (-1, 1),
                            ybounds::Tuple{Real, Real} = (-1, 1),
                            zbounds::Tuple{Real, Real} = (-1, 1))

    # Create points from given intervals
    points = [Point3f0(i, j, k) for i in LinRange(xbounds[1], xbounds[2], 10)
                                for j in LinRange(ybounds[1], ybounds[2], 10)
                                for k in LinRange(zbounds[1], zbounds[2], 10)]

    # Define functions to calculate gradient
    vectorf(v::Vector) = f(v...)
    ∇(p::Point3f0) = ForwardDiff.gradient(vectorf, [p...])

    # Create vectors with specified anchor points
    sumBounds = sum([bounds[2]-bounds[1] for bounds in [xbounds, ybounds, zbounds]])
    vectors = @. ColorVector3D(∇(points), points, 30/sumBounds)
    vectors = differVectors(vectors)

    # Create a plot
    scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], [zbounds[1], zbounds[2]], visible = :false)
    plotVector.(vectors, scene, xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)

    # Display and return the plot
    scene |> display
    return scene

end


"""
    gradientField(f::Function; kwargs...)

Plot a vector field of given two or three-variable function's gradient.

# Arguments
- `showContour::Bool = :false`: flag to indicate whether create contour plot.
- `xbounds::Tuple{Real, Real} = (-1, 1)`: boundaries for x coordinates.
- `ybounds::Tuple{Real, Real} = (-1, 1)`: boundaries for y coordinates.
- `zbounds::Tuple{Real, Real} = (-1, 1)`: boundaries for z coordinates.
"""
function gradientField(f::Function; showContour::Bool = :false,
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        zbounds::Tuple{Real, Real} = (-1, 1))

    # Indicate whether function is two or three-variable
    for method in methods(f)
        if length(method.sig.parameters) == 4
            return gradientField3D(f, xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)
        else
            return gradientField2D(f, xbounds = xbounds, ybounds = ybounds, showContour = showContour)
        end
    end

end


"""
    divergence(point::Array{T, 1}, Fx::Function, Fy::Function [, Fz]) where T <: Real

Compute the divergence in given `point` of the two or three-variable functions by the formula:
    div(F) = ∇ · F
"""
function divergence(point::Array{T, 1}, Fx::Function, Fy::Function,
                    Fz=missing) where T <: Real

    vectorFunctions = []  # functions that take vectors container

    # Define functions for vectors to compute gradient
    vectorFuncX(v::Vector) = Fx(v...)
    vectorFuncY(v::Vector) = Fy(v...)
    push!(vectorFunctions, vectorFuncX)
    push!(vectorFunctions, vectorFuncY)
    if typeof(Fz) <: Function
        vectorFuncZ(v::Vector) = Fz(v...)
        push!(vectorFunctions, vectorFuncZ)
    end

    return sum([ForwardDiff.gradient(func, point)[variable]
                for (variable, func) in enumerate(vectorFunctions)])

end


"""
    curl(point::Array{T, 1}, Fx::Function, Fy::Function [, Fz]) where T <: Real

Compute the curl in given `point` of the two or three-variable functions by the formula:
    curl(F) = ∇ ✖ F
"""
function curl(point::Array{T, 1}, Fx::Function, Fy::Function,
                Fz=missing) where T <: Real

    # Define functions for vectors to compute gradient
    vectorFuncX(v::Vector) = Fx(v...)
    vectorFuncY(v::Vector) = Fy(v...)
    if typeof(Fz) <: Function
        vectorFuncZ(v::Vector) = Fz(v...)
        return [ForwardDiff.gradient(vectorFuncZ, point)[2] - ForwardDiff.gradient(vectorFuncY, point)[3],
                ForwardDiff.gradient(vectorFuncX, point)[3] - ForwardDiff.gradient(vectorFuncZ, point)[1],
                ForwardDiff.gradient(vectorFuncY, point)[1] - ForwardDiff.gradient(vectorFuncX, point)[2]]
    else
        return [0, 0, ForwardDiff.gradient(vectorFuncY, point)[1] - ForwardDiff.gradient(vectorFuncX, point)[2]]
    end

end


# Clear the plot
function clearTrail(body::Body, scene::Scene)

    delete!(scene, scene[end])
    delete!(scene, scene[end])

end


# Add plots for start animation
function dummyPlots2D(body::Body, scene::Scene, xbounds, ybounds)

    lines!(scene, xbounds[1]:xbounds[2], xbounds[1]:xbounds[2]) # dummy trajectory
    lines!(scene, ybounds[1]:ybounds[2], ybounds[1]:ybounds[2]) # dummy dot

end


function dummyPlots3D(body::Body, scene::Scene, xbounds, ybounds)

    lines!(scene, xbounds[1]:xbounds[2], xbounds[1]:xbounds[2], xbounds[1]:xbounds[2]) # dummy trajectory
    lines!(scene, ybounds[1]:ybounds[2], ybounds[1]:ybounds[2], ybounds[1]:ybounds[2]) # dummy dot

end


# Add single frame to animation
function addPlot!(bodies, scene, linewidth; stopFrame = 1)

    clearTrail.(bodies, scene)
    plotBody.(bodies, scene, linewidth, stopFrame = stopFrame)
    return scene

end


"""
    animate(Fx::Function, Fy::Function, startPoints::Array{Array{T, 1}, 1}, filename::String; kwargs...) where T <: Real

Create and save in `filename` an animation of objects' with `startPoints` movement in vector field of given two-variable functions.

# Arguments
- `showField::Bool = :false`: flag to indicate whether plot vector field (arrows).
- `time::Tuple{Real, Real} = (0, 1)`: boundaries for time which positions will be computed in.
- `timePoints::Integer = 100`: number of positions (accuracy) to compute.
- `xbounds::Tuple{Real, Real} = (-1, 1)`: boundaries for x coordinates.
- `ybounds::Tuple{Real, Real} = (-1, 1)`: boundaries for y coordinates.
- `fps::Integer = 30`: framerate of the animation.
- `linewidth::Real = 3`: width of trajectories' lines.
"""
function animate(Fx::Function, Fy::Function,
                    startPoints::Array{Array{T, 1}, 1},
                    filename::String;
                    showField::Bool = :false,
                    time::Tuple{Real, Real} = (0, 1),
                    timePoints::Integer = 100,
                    xbounds::Tuple{Real, Real} = (-1, 1),
                    ybounds::Tuple{Real, Real} = (-1, 1),
                    fps::Integer = 30,
                    linewidth::Real = 3) where T <: Real

    # Create objects
    bodies = @. Body2D(positions(Fx, Fy, startPoints, time = time, timePoints = timePoints))

    # Create plot (with or without vector field)
    if showField
        scene = field(Fx, Fy, xbounds = xbounds, ybounds = ybounds)
    else
        scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], visible = :false)
        for body in bodies
            dummyPlots2D(body, scene, xbounds, ybounds)
        end
    end

    # Make the animation
    record(scene, filename, 1:timePoints, framerate = fps) do frame
        addPlot!(bodies, scene, linewidth, stopFrame = frame)
    end

    return scene

end


"""
    animate(Fx::Function, Fy::Function, Fz::Function, startPoints::Array{Array{T, 1}, 1}, filename::String; kwargs...) where T <: Real

Also for three-variable functions.

# Arguments
- `zbounds::Tuple{Real, Real} = (-1, 1)`: boundaries for z coordinates.
"""
function animate(Fx::Function, Fy::Function, Fz::Function,
                    startPoints::Array{Array{T, 1}, 1}, # starting positions
                    title::String;
                    showField::Bool = :false,
                    time::Tuple{Real, Real} = (0, 1), # time
                    timePoints::Integer = 100,
                    xbounds::Tuple{Real, Real} = (-1, 1),
                    ybounds::Tuple{Real, Real} = (-1, 1),
                    zbounds::Tuple{Real, Real} = (-1, 1),
                    fps::Integer = 30,
                    linewidth::Real = 3) where T <: Real

    # Create objects
    bodies = @. Body3D(positions(Fx, Fy, Fz, startPoints, time = time, timePoints = timePoints))

    # Create plot (with or without vector field)
    if showField
        scene = field(Fx, Fy, Fz, xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)
    else
        scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], [zbounds[1], zbounds[2]], visible = :false)
        for body in bodies
            dummyPlots3D(body, scene, xbounds, ybounds)
        end
    end

    # Make the animation
    record(scene, title, 1:timePoints, framerate = fps) do frame
        addPlot!(bodies, scene, linewidth, stopFrame = frame)
        rotate_cam!(scene, 6/timePoints, 0.0, 0.0) # rotation of a camera
    end

    return scene

end


end
