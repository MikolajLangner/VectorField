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


abstract type Body end


mutable struct Body2D <: Body

    X::Array
    Y::Array
    colors::LinRange
    inBounds::Bool

    Body2D(t::Tuple) = new(t[1], t[2], LinRange(0, length(t[1]), 10length(t[1])), :true)

end


mutable struct Body3D <: Body

    X::Array
    Y::Array
    Z::Array
    colors::LinRange
    inBounds::Bool

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


function cutPositions(body::Body;
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        zbounds::Tuple{Real, Real} = (-1, 1))

    body.X = body.X[xbounds[1] .< body.X .< xbounds[2]]
    body.Y = body.Y[ybounds[1] .< body.Y .< ybounds[2]]
    if typeof(body) == Body3D
        body.Z = body.Z[zbounds[1] .< body.Z .< zbounds[2]]
        # X and Y must have the same length for plotting
        cut = min(length(body.X), length(body.Y), length(body.Z))
        if cut > 0
            body.X = body.X[1:cut]
            body.Y = body.Y[1:cut]
            body.Z = body.Z[1:cut]
            body.colors = LinRange(0, cut, 10*cut)
        else
            body.inBounds = :false
        end
    else
        cut = min(length(body.X), length(body.Y))
        if cut > 0
            body.X = body.X[1:cut]
            body.Y =  body.Y[1:cut]
            body.colors = LinRange(0, cut, 10*cut)
        else
            body.inBounds = :false
        end
    end

    return body

end


function plotBody(body::Body, scene::Scene, linewidth::Real; stopFrame::Integer = length(body.X))

    if body.inBounds
        if typeof(body) == Body3D
            lines!(scene, body.X[1:stopFrame], body.Y[1:stopFrame], body.Z[1:stopFrame],
                    linewidth = linewidth, color = body.colors[1:10stopFrame], colormap = :grayC)
            scatter!(scene, [body.X[stopFrame]], [body.Y[stopFrame]], [body.Z[stopFrame]], markersize=body.X[stopFrame]/20)
        else
            lines!(scene, body.X[1:stopFrame], body.Y[1:stopFrame],
                    linewidth = linewidth, color = body.colors[1:10stopFrame], colormap = :grayC)
            scatter!(scene, [body.X[stopFrame]], [body.Y[stopFrame]], markersize=linewidth/70)
        end
    end

    return scene

end


# function for calculating position (x(t), y(t)) of an object
function positions2D(Fx::Function, Fy::Function, startPoint::Array{T, 1} = [0, 0]; # starting position
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


# function for plotting trajectory based on calculations from "position2D"
function trajectory2D(Fx::Function, Fy::Function,
                        startPoints::Array{Array{T, 1}, 1};  # starting positions
                        showField::Bool = :false,
                        time::Tuple{Real, Real} = (0, 1), # time
                        timePoints::Integer = 100,
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        linewidth::Real = 3) where T <: Real

    bodies = @. Body2D(positions2D(Fx, Fy, startPoints, time = time, timePoints = timePoints))
    bodies = cutPositions.(bodies, xbounds = xbounds, ybounds = ybounds)

    if showField
        scene = field2D(Fx, Fy, xbounds = xbounds, ybounds = ybounds)
    else
        scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], visible = :false)
    end

    plotBody.(bodies, scene, linewidth)

    scene |> display
    return scene

end


# function for calculating position (x(t), y(t), z(t)) of an object
function positions3D(Fx::Function, Fy::Function, Fz::Function,
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

    prob = ODEProblem(diffeqs, Float64.(startPoint), Float64.(time))
    sol = solve(prob)
    ts = LinRange(time[1], time[2], timePoints)
    X, Y, Z = sol(ts, idxs = 1), sol(ts, idxs = 2), sol(ts, idxs = 3)
    return X, Y, Z # X - Array of x coords, Y, Z similarly

end


# function for plotting trajectory based on calculations from "position3D"
function trajectory3D(Fx::Function, Fy::Function, Fz::Function,
                        startPoints::Array{Array{T, 1}, 1};
                        showField::Bool = :false,
                        time::Tuple{Real, Real} = (0, 1),
                        timePoints::Integer = 100,
                        xbounds::Tuple{Real, Real} = (-1, 1),
                        ybounds::Tuple{Real, Real} = (-1, 1),
                        zbounds::Tuple{Real, Real} = (-1, 1),
                        linewidth::Real = 3) where T <: Real

    bodies = @. Body3D(positions3D(Fx, Fy, Fz, startPoints, time = time, timePoints = timePoints))
    bodies = cutPositions.(bodies, xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)

    if showField
        scene = field3D(Fx, Fy, Fz,
                    xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)
    else
        scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], [zbounds[1], zbounds[2]] visible = :false)
    end

    plotBody.(bodies, scene, linewidth)
    scene |> display
    return scene

end


# Plot a vector field from gradient of given two-variables function
function gradientField2D(f::Function,
                            xbounds::Tuple{Real, Real},
                            ybounds::Tuple{Real, Real},
                            showContour::Bool = :false)

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

end


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


function clearTrail(body::Body, scene::Scene)
    delete!(scene, scene[end])
    delete!(scene, scene[end])
end

function dummyPlots2D(body::Body, scene::Scene, xbounds, ybounds)
    lines!(scene, xbounds[1]:xbounds[2], xbounds[1]:xbounds[2]) # dummy trajectory
    lines!(scene, ybounds[1]:ybounds[2], ybounds[1]:ybounds[2]) # dummy dot
end

function dummyPlots3D(body::Body, scene::Scene, xbounds, ybounds)
    lines!(scene, xbounds[1]:xbounds[2], xbounds[1]:xbounds[2], xbounds[1]:xbounds[2]) # dummy trajectory
    lines!(scene, ybounds[1]:ybounds[2], ybounds[1]:ybounds[2], ybounds[1]:ybounds[2]) # dummy dot
end

function addPlot!(bodies, scene, linewidth; stopFrame = 1)

    clearTrail.(bodies, scene)
    plotBody.(bodies, scene, linewidth, stopFrame = stopFrame)
    return scene

end


function animate2D(Fx::Function, Fy::Function,
                    startPoints::Array{Array{T, 1}, 1}, # starting positions
                    title::String;
                    showField::Bool = :false,
                    time::Tuple{Real, Real} = (0, 1), # time
                    timePoints::Integer = 100,
                    xbounds::Tuple{Real, Real} = (-1, 1),
                    ybounds::Tuple{Real, Real} = (-1, 1),
                    fps::Integer = 30,
                    linewidth::Real = 3) where T <: Real

    bodies = @. Body2D(positions2D(Fx, Fy, startPoints, time = time, timePoints = timePoints))

    if showField
        scene = field2D(Fx, Fy, xbounds = xbounds, ybounds = ybounds)
    else
        scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], visible = :false)
        for body in bodies
            dummyPlots2D(body, scene, xbounds, ybounds)
        end
    end

    record(scene, title, 1:timePoints, framerate = fps) do frame
        addPlot!(bodies, scene, linewidth, stopFrame = frame)
    end

    return scene

end


function animate3D(Fx::Function, Fy::Function, Fz::Function,
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

    bodies = @. Body3D(positions3D(Fx, Fy, Fz, startPoints, time = time, timePoints = timePoints))

    if showField
        scene = field3D(Fx, Fy, Fz, xbounds = xbounds, ybounds = ybounds, zbounds = zbounds)
    else
        scene = lines([xbounds[1], xbounds[2]], [ybounds[1], ybounds[2]], [zbounds[1], zbounds[2]], visible = :false)
        for body in bodies
            dummyPlots3D(body, scene, xbounds, ybounds)
        end
    end

    record(scene, title, 1:timePoints, framerate = fps) do frame
        addPlot!(bodies, scene, linewidth, stopFrame = frame)
        rotate_cam!(scene, 6/timePoints, 0.0, 0.0)
    end

    return scene

end


# end
