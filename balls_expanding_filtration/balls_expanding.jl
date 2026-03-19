using Plots
using LinearAlgebra
using Ripserer   # すでに読み込んでいる想定


# 円盤
function circle(x, y, r; n=60)
    θ = range(0, 2π; length=n)
    Shape(x .+ r*cos.(θ), y .+ r*sin.(θ))
end
dist(p, q) = hypot(p[1]-q[1], p[2]-q[2])

# Čech 辺（1-skeleton）
function cech_edges(points, r)
    n = length(points)
    edges = Tuple{Int,Int}[]
    thr = 2r
    for i in 1:n-1
        pi = points[i]
        for j in i+1:n
            if dist(pi, points[j]) <= thr
                push!(edges, (i, j))
            end
        end
    end
    return edges
end

# 三角形（2-simplex）判定：外接円半径 <= r
function circumradius(p, q, s)
    a = dist(q, s)
    b = dist(p, s)
    c = dist(p, q)
    area2 = abs((q[1]-p[1])*(s[2]-p[2]) - (q[2]-p[2])*(s[1]-p[1]))
    area2 == 0.0 && return Inf
    area = area2 / 2
    return (a*b*c) / (4*area)
end

function cech_triangles(points, r; max_triangles=typemax(Int))
    n = length(points)
    tris = NTuple{3,Int}[]
    thr = 2r

    adj = [BitSet() for _ in 1:n]
    for i in 1:n-1
        for j in i+1:n
            if dist(points[i], points[j]) <= thr
                push!(adj[i], j); push!(adj[j], i)
            end
        end
    end

    for i in 1:n-2
        for j in adj[i]
            j <= i && continue
            for k in intersect(adj[i], adj[j])
                k <= j && continue
                if circumradius(points[i], points[j], points[k]) <= r
                    push!(tris, (i, j, k))
                    length(tris) >= max_triangles && return tris
                end
            end
        end
    end
    return tris
end


# VR(=Rips) 2-単体（三角形）
# t: Rips の距離しきい値（dist <= t なら辺がある）
function rips_triangles(points, t; max_triangles=typemax(Int))
    n = length(points)
    tris = NTuple{3,Int}[]

    # まず隣接（辺があるか）を作って枝刈り
    adj = [BitSet() for _ in 1:n]
    for i in 1:n-1
        pi = points[i]
        for j in i+1:n
            if dist(pi, points[j]) <= t
                push!(adj[i], j)
                push!(adj[j], i)
            end
        end
    end

    # 三角形 (i,j,k): j,k が i の近傍で、かつ k が j の近傍
    for i in 1:n-2
        for j in adj[i]
            j <= i && continue
            for k in intersect(adj[i], adj[j])
                k <= j && continue
                push!(tris, (i, j, k))
                if length(tris) >= max_triangles
                    return tris
                end
            end
        end
    end

    return tris
end

# Ripserer の diagram から (birth, death)
function intervals_from(diagram)
    bs = Float64.(getproperty.(diagram, :birth))
    ds = Float64.(getproperty.(diagram, :death))
    return collect(zip(bs, ds))
end

# “伸びるバーコード”
function plot_growing_barcode!(p, intervals, t; y0=0, dy=1.0, lc=:black, lw=3)
    y = y0
    for (b, d) in intervals
        b > t && continue
        right = isfinite(d) ? min(t, d) : t
        right <= b && continue
        plot!(p, [b, right], [y, y]; lc=lc, lw=lw, label=false)
        y += dy
    end
end

# 上位k本を選ぶ（長さ順）
# death=Inf は長さ=Inf 扱い。death が tmax より大きい場合は tmax までで評価。
function select_topk_intervals(intervals, k::Union{Int,Nothing}, tmax::Float64)
    k === nothing && return intervals
    k <= 0 && return Tuple{Float64,Float64}[]

    lens = map(intervals) do (b, d)
        if !isfinite(d)
            Inf
        else
            max(0.0, min(d, tmax) - b)
        end
    end
    idx = sortperm(lens; rev=true)
    k = min(k, length(idx))
    return intervals[idx[1:k]]
end




function make_filtration_plus_barcode_gif(points, result;
        rmax::Float64,
        nframes::Int=80,
        fps::Int=10,
        outpath::AbstractString="filtration_barcode.gif",
        draw_triangles::Bool=false,
        max_triangles::Int=2000,
        ball_alpha::Float64=0.25,
        tri_alpha::Float64=0.15,
        topk0::Union{Int,Nothing}=nothing,
        topk1::Union{Int,Nothing}=nothing
    )

    x = first.(points); y = last.(points)

    # --- バーコード用区間（r軸に変換：t/2）---
    ints0_all = intervals_from(result[1])
    ints1_all = intervals_from(result[2])
    ints0_all = [(b/2, d/2) for (b, d) in ints0_all]
    ints1_all = [(b/2, d/2) for (b, d) in ints1_all]

    # --- PD表示用（r軸で表示したいので、こちらも t/2 に変換した“図”を作る）---
    # ここは「点の散布図」で描くのが確実（Ripsererのplot結果を直接いじらない）
    dgm0 = ints0_all
    dgm1 = ints1_all

    # 表示範囲（データ）
    xmin, xmax = minimum(x), maximum(x)
    ymin, ymax = minimum(y), maximum(y)
    pad = 0.1 * max(xmax - xmin, ymax - ymin)
    xlims = (xmin - pad - rmax, xmax + pad + rmax)
    ylims = (ymin - pad - rmax, ymax + pad + rmax)

    # バーコード・PD横軸は r
    tmax = rmax

    # 間引き
    ints0 = select_topk_intervals(ints0_all, topk0, tmax)
    ints1 = select_topk_intervals(ints1_all, topk1, tmax)

    n0 = length(ints0)
    n1 = length(ints1)

    title_suffix = ""
    if topk0 !== nothing || topk1 !== nothing
        s0 = (topk0 === nothing) ? "H0:All" : "H0:Top$(topk0)"
        s1 = (topk1 === nothing) ? "H1:All" : "H1:Top$(topk1)"
        title_suffix = "($s0, $s1)"
    end

    # PDの表示範囲（r軸なので [0, rmax] あたりに合わせる）
    pd_xlim = (0, tmax)
    pd_ylim = (0, tmax)


    # --- 右上：静的な Persistence diagram（r軸） ---
    p_pd_static = plot(
        xlims=(0, tmax), ylims=(0, tmax),
        aspect_ratio=1,
        xlabel="birth (r)", ylabel="death (r)",
        legend=false,
        title="Persistence diagram (H0 black, H1 red)"
    )

scatter!(p_pd_static, first.(dgm0), last.(dgm0); ms=3,msw=0, c=:black)
scatter!(p_pd_static, first.(dgm1), last.(dgm1); ms=3,msw=0, c=:red)
plot!(p_pd_static, [0, tmax], [0, tmax]; lc=:gray, lw=1, alpha=0.6, label=false)

    anim = @animate for r in range(0, rmax; length=nframes)
        t = r

        # --- 左：ボール膨張 ---
        p1 = plot(aspect_ratio=1, xlims=xlims, ylims=ylims, legend=false)
        scatter!(p1, x, y; ms=2)

        for (xi, yi) in zip(x, y)
            plot!(p1, circle(xi, yi, r); fc=:blue, fa=ball_alpha, lc=:blue, lw=0, label=false)
        end

        for (i, j) in cech_edges(points, r)
            plot!(p1, [x[i], x[j]], [y[i], y[j]]; lc=:black, lw=1, alpha=0.8, label=false)
        end

        if draw_triangles
            tris = rips_triangles(points, 2r; max_triangles=max_triangles)
            for (i, j, k) in tris
                poly = Shape([x[i], x[j], x[k]], [y[i], y[j], y[k]])
                plot!(p1, poly; lw=0, lc=:red, fillcolor=:red, fillalpha=tri_alpha, label=false)
            end
        end

        # 表示位置
        annotate!(p1,
            xlims[1] - 0.25*(xlims[2]-xlims[1]),
            ylims[2] - 0*(ylims[2]-ylims[1]),
            text("r = $(round(r, digits=4))", 15)
        )

        p_pd = deepcopy(p_pd_static)   # 安全のため（毎フレーム描画状態が汚れない）

        # H0: black, H1: red
        scatter!(p_pd, first.(dgm0), last.(dgm0); ms=2, c=:black)
        scatter!(p_pd, first.(dgm1), last.(dgm1); ms=2, c=:red)

        # 対角線
        plot!(p_pd, [0, tmax], [0, tmax]; lc=:gray, lw=1, alpha=0.6, label=false)

        # --- 下：バーコード ---
        gap = 3
        offset = n0 + gap
        p2 = plot(xlims=(0, tmax), ylims=(-1, n0 + n1 + gap + 2),
                  legend=false, yticks=false, xlabel="r",
                  title="Barcodes (H0:black, H1:red), Rips $title_suffix")
        p2 = plot(
    xlims=(0, tmax), ylims=(-1, n0 + n1 + gap + 2),
    legend=false, yticks=false, xlabel="r",
    title="Barcodes  (H0:black, H1:red), Rips $title_suffix",
    left_margin=12Plots.mm,
    right_margin=12Plots.mm
)

        plot_growing_barcode!(p2, ints0, t; y0=0,      dy=1.0, lc=:black, lw=3)
        plot_growing_barcode!(p2, ints1, t; y0=offset, dy=1.0, lc=:red,   lw=3)

        # --- 合成：上段(左=データ, 右=PD)、下段(バーコードを横幅いっぱい) ---
        plot(p1,  p_pd, p2;
     layout = @layout([a b; c ]),
     size   = (1200, 700))
    end

    gif(anim, outpath; fps=fps)
    return outpath
end

