module VectorField


using Makie
using LinearAlgebra


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
    for vector in vectors
        arrows!([vector.anchorPoint], [vector.direction],
                arrowsize=0.05vector.size, linewidth=2vector.size,
                linecolor=vector.color, arrowcolor=vector.color)
    end

    # Display and return the plot
    scene |> display
    return scene

end


# Example
field2D((x,y)->[x^2+y^2,x-y^2],(-1,1),(-1,1))


# Main fuction for 3D plotting
function field3D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real}, zbounds::Tuple{Real,Real})

    # Create points in ℝ³
    points = [Point3f0(i, j, k) for i in LinRange(xbounds[1], xbounds[2], 10)
                                for j in LinRange(ybounds[1], ybounds[2], 10)
                                for k in LinRange(zbounds[1], zbounds[2], 10)]

    # Define new function; `f` that takes points in ℝ³
    pointf(p::Point3f0) = f(p[1], p[2], p[3])

    # Create vectors
    vectors = @. ColorVector3D(pointf(points), points, 10)
    vectors = differVectors(vectors)

    # Create a plot
    scene = Scene()
    for vector in vectors
        arrows!([vector.anchorPoint], [vector.direction],
                arrowsize=0.05vector.size, linewidth=2vector.size,
                linecolor=vector.color, arrowcolor=vector.color)
    end

    # Display and return the plot
    scene |> display
    return scene

end


# Example
field3D((x,y,z)->[x^2+y^2,x-y,z],(-3,2),(-1,2),(-1,2))


end
