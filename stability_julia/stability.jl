using Ripserer
using Plots
using Random
MT = MersenneTwister(12345)
MT2 = MersenneTwister(123456)
rand(MT)

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

points = noisy_circle(100; noise=0.2)
pointb = noisy_circle(100; noise=0.2)

for i in 1:100
    push!(points, (pointb[i][1] +1.2, pointb[i][2]))
end

##########animation
rng = MersenneTwister(1234)
function add_uniform_noise(points::Vector{Tuple{Float64,Float64}}; ε::Float64=0.03, rng=Random.default_rng())
    [(x + (2rand(rng)-1)*ε, y + (2rand(rng)-1)*ε) for (x,y) in points]
end

for i in 1:length(points)
        noise = add_uniform_noise(points)
        result = ripserer(noise)[2]
        plt_pts = scatter(
        points;
        legend=false,
        aspect_ratio=1,
        xlim=(-2.2, 2.2),
        ylim=(-2.2, 2.2),
        title="Data",
    )
    plt_diag = plot(result; infinity=3)

    plot(plt_pts, plt_diag; size=(800, 400))
end
noise = add_uniform_noise(points)

xt_pts = -2:1:2
yt_pts = -2:1:2

# PDの軸を固定（infinity=3 に合わせて 0〜3）
pd_min, pd_max = 0.0, 1.5
xt_pd = 0:0.5:3
yt_pd = 0:0.5:3

anim = @animate for i in 1:length(points)
    noisy_points = add_uniform_noise(points; ε=0.05, rng=rng)
    result = ripserer(noisy_points)[2]
    plt_pts =scatter(noisy_points;  title="Data with noise", markersize=2, aspect_ratio=1,label = false)
 
plt_pts = scatter(
        first.(noisy_points), last.(noisy_points);
        legend=false,
        aspect_ratio=1,
        xlim=(-1.2, 2.5),
        ylim=(-1.2, 1.5),
        xticks=xt_pts,
        yticks=yt_pts,
        title="Data with noise",
         markersize=2,
    )

    plt_diag = plot(
        result;
        infinity=3,
        xlim=(pd_min, pd_max),
        ylim=(pd_min, pd_max),
        xticks=xt_pd,
        yticks=yt_pd,
        aspect_ratio=1,
        title="Persistence Diagram",
    )

    plot(plt_pts, plt_diag; size=(800, 400))
end

gif(anim, "noisy_pd.gif"; fps=10)



