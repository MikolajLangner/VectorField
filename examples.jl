include("vector_field.jl")
using .VectorField

# 3D fields

field((x, y, z)->-y, (x, y, z)->x, (x, y, z)->sin(z))

field((x, y, z)->x-y-z, (x, y, z)->y-x-z, (x, y, z)->z-x-y,
        xbounds = (-2, 3), ybounds = (-1, 4), zbounds = (0, 2.5))

gradientField((x, y, z)->x+sin(y)-cos(z), ybounds = (0, 2pi), zbounds = (-pi, pi))


# animations

animate((x, y)->sin(exp(-y)), (x, y)->cos(exp(x)), [[i, (-1)^j * sqrt(1-i^2)] for i in -1:0.2:1 for j in 0:1], "trigexp.gif",
                showField = :true, timePoints = 60, xbounds = (-pi, pi), ybounds = (-pi, pi), time = (0, 2))

animate((x, y, z)->cos(1/exp(z)), (x, y, z)->sin(1/exp(z)), (x, y, z)->1/(1+x^2+y^2),
        [[i, (-1)^j * sqrt(1-i^2), k] for i in -1:0.5:1 for j in 0:1 for k in -1:0.5:1], "3dtrigexp.gif",
        xbounds = (-3, 3), ybounds = (-3, 3), zbounds = (-1, 4))
