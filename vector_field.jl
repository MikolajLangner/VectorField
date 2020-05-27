module VectorField


using Makie
using LinearAlgebra


function field2D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real})
    points = [Point2f0(i, j) for i in LinRange(xbounds[1], xbounds[2], 20)
                             for j in LinRange(ybounds[1], ybounds[2], 20)]
    pointf(p::Point2f0) = f(p[1], p[2])
    vectors = @. Point2f0(normalize(pointf(points))/20)
    scene = arrows(points, vectors, arrowsize=0.02, linewidth=1, linecolor=:red)
    scene |> display
    return scene
end


function field3D(f::Function, xbounds::Tuple{Real, Real}, ybounds::Tuple{Real, Real}, zbounds::Tuple{Real,Real})
    points = [Point3f0(i, j, k) for i in LinRange(xbounds[1], xbounds[2], 10)
                             for j in LinRange(ybounds[1], ybounds[2], 10)
                             for k in LinRange(zbounds[1], zbounds[2], 10)]
    pointf(p::Point3f0) = f(p[1], p[2], p[3])
    vectors = @. Point3f0(normalize(pointf(points))/10)
    scene = arrows(points, vectors, arrowsize=0.02, linewidth=1, linecolor=:black)
    scene |> display
    return scene
end


end
