module VectorField


using Makie
using LinearAlgebra


abstract type ColorVector end


mutable struct ColorVector2D <: ColorVector

    magnitude::Float64
    direction::Point2f0
    anchorPoint::Point2f0
    color::RGBf0
    size::Float64

    ColorVector2D(vector::Array{T, 1}, point::Point2f0, decreaser::Real) where T <: Real = new(norm(vector),
                                                                        Point2f0(normalize(vector)/decreaser),
                                                                        point, RGBf0(0, 0, 0), 1.0)

end


mutable struct ColorVector3D <: ColorVector

    magnitude::Float64
    direction::Point3f0
    anchorPoint::Point3f0
    color::RGBf0
    size::Float64

    ColorVector3D(vector::Array{T, 1}, point::Point3f0, decreaser::Real) where T <: Real = new(norm(vector),
                                                                        Point3f0(normalize(vector)/decreaser),
                                                                        point, RGBf0(0, 0, 0), 1.0)

end


Base.isless(v::ColorVector, w::ColorVector) = v.magnitude < w.magnitude ? true : false
Base.isequal(v::ColorVector, w::ColorVector) = v.magnitude == w.magnitude ? true : false


function changeLength(v::ColorVector, longest::Real)

    v.size = v.magnitude / longest
    v.direction = v.direction * v.size
    return v

end


function changeColor(v::ColorVector, longest::Real)

    v.color = v.direction[end] >= 0 ? RGBf0(v.size, 1-v.size, 0) : RGBf0(0, 1-v.size, v.size)
    return v

end


function differVectors(vectors::Array{T, 1})::Array{T, 1} where T <: ColorVector

    vectors = changeLength.(vectors, vectors[end].magnitude)
    vectors = changeColor.(vectors, vectors[end].magnitude)
    return vectors

end


function field2D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real})

    points = [Point2f0(i, j) for i in LinRange(xbounds[1], xbounds[2], 20)
                             for j in LinRange(ybounds[1], ybounds[2], 20)]

    pointf(p::Point2f0) = f(p[1], p[2])
    vectors = @. ColorVector2D(pointf(points), points, 10)
    vectors = differVectors(vectors)

    scene = Scene()
    for vector in vectors
        arrows!([vector.anchorPoint], [vector.direction],
                arrowsize=0.05vector.size, linewidth=2vector.size,
                linecolor=vector.color, arrowcolor=vector.color)
    end

    scene |> display
    return scene

end


# Example
field2D((x,y)->[x^2+y^2,x-y^2],(-1,1),(-1,1))


function field3D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real}, zbounds::Tuple{Real,Real})

    points = [Point3f0(i, j, k) for i in LinRange(xbounds[1], xbounds[2], 10)
                                for j in LinRange(ybounds[1], ybounds[2], 10)
                                for k in LinRange(zbounds[1], zbounds[2], 10)]

    pointf(p::Point3f0) = f(p[1], p[2], p[3])
    vectors = @. ColorVector3D(pointf(points), points, 10)
    vectors = differVectors(vectors)

    scene = Scene()
    for vector in vectors
        arrows!([vector.anchorPoint], [vector.direction],
                arrowsize=0.05vector.size, linewidth=2vector.size,
                linecolor=vector.color, arrowcolor=vector.color)
    end

    scene |> display
    return scene

end


# Example
field3D((x,y,z)->[x^2+y^2,x-y,z],(-3,2),(-1,2),(-1,2))


end
