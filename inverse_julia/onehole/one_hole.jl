using LinearAlgebra
using Plots
using Ripserer
using Random

Random.seed!(1337)

# -----------------------------
# 1. データ生成（Annulus）
# -----------------------------
function annulus(n; r1=1.0, r2=2.0, offset=(0.0, 0.0))
    pts = Tuple{Float64,Float64}[]
    while length(pts) < n
        p = 2r2 .* rand(2) .- r2
        if r1 < norm(p) < r2
            push!(pts, (p[1] + offset[1], p[2] + offset[2]))
        end
    end
    return pts
end

data = annulus(300)

scatter(data; label="data", ms=2, aspect_ratio=1)
savefig("data_without_cycle.png")

# -----------------------------
# 2. パーシステンス図（H1）
# -----------------------------
res_pd = ripserer(data)          # reps=false (default)
dgm1   = res_pd[2]               # H1 diagram (Ripserer は 1-indexed: [1]=H0, [2]=H1)

plot(dgm1)
savefig("1PD.png")

# -----------------------------
# 3. 逆解析：代表（co)cycle を使って「穴」を再構成
#    - reps=true で cocycle representative を取得
# -----------------------------
res_rep = ripserer(data; reps=true)
dgm1_rep = res_rep[2]
filt_rep = dgm1_rep.filtration

# 「対角線から遠い点」= よく残る（長い）区間を選ぶ：最後が最長になっていることが多い
most_persistent = dgm1_rep[end]

# 例：2つ目の候補（インデックスで選んでいた部分を、そのまま残す）
# ここはデータや結果でインデックスが変わるので、必要なら “長さ順で2番目” に変えるのがおすすめ
second_candidate = dgm1_rep[48]

# -----------------------------
# 4. 再構成するスケールの選び方
#    - birth 直後 / 区間の中点 など
# -----------------------------
function midpoint_scale(interval)
    # birth/death は Ripserer が提供
    b = birth(interval)
    d = death(interval)
    # death が Inf の場合は適当な扱いが必要（ここではそのまま）
    return (b + d) / 2
end

# 中点スケール（
t_most   = midpoint_scale(most_persistent)
t_second = midpoint_scale(second_candidate)

# -----------------------------
# 5. サイクル再構成 & 可視化
# -----------------------------
cycle_most   = reconstruct_cycle(filt_rep, most_persistent, t_most)
cycle_second = reconstruct_cycle(filt_rep, second_candidate, t_second)

scatter(data; label="data", ms=2, aspect_ratio=1)
plot!(cycle_most, data; label="cycle (most persistent)")
savefig("data_hole_a.png")

scatter(data; label="data", ms=2, aspect_ratio=1)
plot!(cycle_second, data; label="cycle (second candidate)")
savefig("data_second_hole_a.png")