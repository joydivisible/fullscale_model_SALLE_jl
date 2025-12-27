#!/usr/bin/env julia
# Map experimental feeds (mass-basis) to Aspen tie-lines (mass-fraction LLE)
using CSV, DataFrames, LinearAlgebra, Printf

println("=== Tie-line mapping for experimental feeds (mass fraction basis) ===")

# --- LOAD DATA ---
exp_feeds = CSV.read("pre-processed in data.csv", DataFrame)
aspen_tie = CSV.read("../TrainingSets/DATA_SALT_2_5.csv", DataFrame)

# --- COMPUTE MASS FRACTIONS OF EXPERIMENTAL FEEDS ---
mW = exp_feeds[!, "Amount Water (g)"]
mF = exp_feeds[!, "Amount Furfural (g)"]
mEA = exp_feeds[!, "Amount Ethyl Acetate (g)"]

mtot = mW .+ mF .+ mEA
zW = mW ./ mtot
zF = mF ./ mtot
zEA = mEA ./ mtot

exp_feed_fracs = hcat(zW, zF, zEA)

# --- ASPEN FEED MASS FRACTIONS (normalize flows) ---
flowW = aspen_tie[!, "FLOWWA"]
flowF = aspen_tie[!, "FLOWFUR"]
flowEA = aspen_tie[!, "FLOWEA"]

total_flow = flowW .+ flowF .+ flowEA
feed_aspen = hcat(flowW ./ total_flow, flowF ./ total_flow, flowEA ./ total_flow)


# --- FUNCTION: Compute closest or interpolated tie-line ---
function find_tieline_for_feed(zfeed::Vector{Float64}, feed_aspen::Matrix{Float64}, tie_df::DataFrame)
    # distance in ternary space
    dists = [norm(zfeed .- feed_aspen[i, :]) for i in 1:size(feed_aspen, 1)]
    idx_sorted = sortperm(dists)
    i1, i2 = idx_sorted[1:2]
    d1, d2 = dists[i1], dists[i2]

    if d1 < 0.005
        # good match (convert DataFrameRow to Vector)
        wr = collect(values(tie_df[i1, ["XWIN2","XFIN2","XEAIN2"]]))  # water-rich
        ea = collect(values(tie_df[i1, ["XWIN1","XFIN1","XEAIN1"]]))  # EA-rich
        return (wr, ea, i1, d1, 1.0)
    else
        # weighted linear interpolation by inverse distance
        w1 = 1/d1
        w2 = 1/d2
        α = w2 / (w1 + w2)

        wr1 = collect(values(tie_df[i1, ["XWIN2","XFIN2","XEAIN2"]]))
        wr2 = collect(values(tie_df[i2, ["XWIN2","XFIN2","XEAIN2"]]))
        ea1 = collect(values(tie_df[i1, ["XWIN1","XFIN1","XEAIN1"]]))
        ea2 = collect(values(tie_df[i2, ["XWIN1","XFIN1","XEAIN1"]]))

        wr = (1 - α) .* wr1 .+ α .* wr2
        ea = (1 - α) .* ea1 .+ α .* ea2

        return (wr, ea, (i1, i2), d1, α)
    end
end


# --- LOOP OVER EXPERIMENTAL FEEDS ---
mapped = DataFrame(
    Run = Int[],
    Feed_W = Float64[], Feed_F = Float64[], Feed_EA = Float64[],
    WR_W = Float64[], WR_F = Float64[], WR_EA = Float64[],
    EA_W = Float64[], EA_F = Float64[], EA_EA = Float64[],
    TieLineRef = String[], Distance = Float64[]
)

for (i, row) in enumerate(eachrow(exp_feed_fracs))
    zfeed = collect(row)
    wr, ea, ref, dist, α = find_tieline_for_feed(zfeed, feed_aspen, aspen_tie)

    push!(mapped, (
        i,
        zfeed[1], zfeed[2], zfeed[3],
        wr[1], wr[2], wr[3],
        ea[1], ea[2], ea[3],
        string(ref),
        dist
    ))
end

CSV.write("mapped_tielines_aspen_massbasis.csv", mapped)
println("💾 Written mapped_tielines_aspen_massbasis.csv with interpolated endpoints.")

# --- Summary print ---
for r in eachrow(mapped)
    @printf("Run %d: Feed=(%.4f, %.4f, %.4f) → Water-rich=(%.4f, %.4f, %.4f) | EA-rich=(%.4f, %.4f, %.4f) [%s, d=%.4f]\n",
        r.Run, r.Feed_W, r.Feed_F, r.Feed_EA,
        r.WR_W, r.WR_F, r.WR_EA, r.EA_W, r.EA_F, r.EA_EA, r.TieLineRef, r.Distance)
end

println("\n✅ Tie-line mapping complete (mass-basis). Each experimental feed now has Aspen-derived equilibrium endpoints.")
