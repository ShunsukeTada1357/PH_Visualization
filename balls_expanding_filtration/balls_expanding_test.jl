##########################################
using Ripserer
using Plots
using Random
MT = MersenneTwister(12345)
MT2 = MersenneTwister(123456)

function noisy_circle(n; r1=1, r2=1, noise=0.1)
    points = NTuple{2,Float64}[]
    for _ in 1:n
        θ = 2π * rand(MT)
        point = (
            r1 * sin(θ) + noise * rand(MT) - noise / 2,
            r2 * cos(θ) + noise * rand(MT2) - noise / 2,
        )
        push!(points, point)
    end
    return points
end
points = noisy_circle(20; noise=0.2)
result = ripserer(points; reps=false)

make_filtration_plus_barcode_gif(points, result; rmax=1.2,draw_triangles = true, outpath="all.gif")

make_filtration_plus_barcode_gif(points, result; rmax=0.2, topk0=10, topk1=10, outpath="topk.gif")

#######################################