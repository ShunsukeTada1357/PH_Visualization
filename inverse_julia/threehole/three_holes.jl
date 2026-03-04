using Ripserer
using Plots
using Random

# -----------------------------
# 1) データ生成：ノイズ付き楕円（円）点群
# -----------------------------
function noisy_circle(rngθ::AbstractRNG, rngnoise::AbstractRNG, n;
        r1=1.0, r2=1.0, noise=0.1
    )
    pts = NTuple{2,Float64}[]
    for _ in 1:n
        θ = 2π * rand(rngθ)
        push!(pts, (
            r1 * sin(θ) + noise * (rand(rngnoise) - 0.5),
            r2 * cos(θ) + noise * (rand(rngnoise) - 0.5),
        ))
    end
    return pts
end

# 2つの円を合成（片方を平行移動）
function two_noisy_circles(; n=100, noise=0.2, shift=(1.2, 0.0), seed1=12345, seed2=123456)
    rngθ  = MersenneTwister(seed1)
    rngn1 = MersenneTwister(seed1 + 1)
    rngn2 = MersenneTwister(seed2)

    A = noisy_circle(rngθ, rngn1, n; noise=noise)
    B = noisy_circle(rngθ, rngn2, n; noise=noise)
    B = [(p[1] + shift[1], p[2] + shift[2]) for p in B]
    return vcat(A, B)
end

points = two_noisy_circles(; n=100, noise=0.2, shift=(1.2, 0.0))

# -----------------------------
# 2) PD描画（H1）とデータ散布図
# -----------------------------
res_rep = ripserer(points; reps=true)   # representative 付き
dgm1    = res_rep[2]                    # H1
filt    = dgm1.filtration               # reconstruct_cycle に使う filtration

plt_pts = scatter(points;
    legend=false, aspect_ratio=1,
    xlim=(-2.2, 2.2), ylim=(-2.2, 2.2),
    title="Data"
)

plt_pd = plot(dgm1; infinity=3, title="H1 persistence diagram")
plot(plt_pts, plt_pd; size=(900, 420))
savefig("two_cycles_data_and_PD.png")

plot(dgm1; infinity=3)
savefig("two_cycles_PD.png")

# -----------------------------
# 3) 逆解析：H1の「長い順」上位k本を選び、サイクル再構成して保存
# -----------------------------
# 区間長（persistence）が長い順に上位kを取る（death=Infは上位になりがちなので必要なら除外）
function topk_intervals_by_persistence(diagram, k; ignore_infinite=true)
    ints = collect(diagram)
    if ignore_infinite
        ints = [iv for iv in ints if isfinite(death(iv))]
    end
    sort!(ints; by=iv -> (death(iv) - birth(iv)), rev=true)
    return ints[1:min(k, length(ints))]
end

# 区間の「中点スケール」： (birth + death)/2
midpoint_scale(iv) = (birth(iv) + death(iv)) / 2

topk = topk_intervals_by_persistence(dgm1, 3)

for (idx, iv) in enumerate(topk)
    t = midpoint_scale(iv)
    cyc = reconstruct_cycle(filt, iv, t)

    p = scatter(points; label="data", ms=2, aspect_ratio=1,
                title="Cycle $(idx)  (birth=$(round(birth(iv),digits=3)), death=$(round(death(iv),digits=3)))")
    plot!(p, cyc, points; label="cycle_$(idx)")
    savefig(p, "cycle_$(idx).png")
end

# 参考：データだけ
scatter(points; legend=false, ms=2, aspect_ratio=1)
savefig("two_cycles.png")