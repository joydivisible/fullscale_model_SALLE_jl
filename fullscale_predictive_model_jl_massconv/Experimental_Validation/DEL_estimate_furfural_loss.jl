#!/usr/bin/env julia
using CSV, DataFrames, Plots, StatsPlots
include("../calculate_phase_split.jl")


# Load data
ternary_df = CSV.read("ternary_compositions.csv", DataFrame)
feed_df = CSV.read("pre-processed out data.csv", DataFrame)
τ_fitted = read_tau_matrix("../fitted_tau_matrix_2_5.csv")

# Simulation and results
results = DataFrame(Run = Int[], Exp_FUR_mgL = Float64[], Sim_FUR_mgL = Float64[])
rho_org = 902.0 # g/L (EA density)

for run_num in ternary_df.Run
    row = ternary_df[ternary_df.Run .== run_num, :]
    feed_row = feed_df[feed_df."Run #" .== run_num, :]
    
    # Convert percentages to mass fractions
    z_feed = [
        feed_row."Water Feed, mass %"[1] / 100,
        feed_row."Furfural Feed, mass %"[1] / 100,
        feed_row."Ethyl Acetate Feed, mass %"[1] / 100
    ]
    
    # Get experimental furfural in EA-rich phase (already in mg/L)
    fur_in_EA_str = split(feed_row."Furfural in EA-rich"[1], "/")[1]
    exp_furfural_mgL = parse(Float64, fur_in_EA_str)

    # Solve LLE
    try
        phase1, phase2 = solve_lle(z_feed, τ_fitted)
        model_EA_rich = phase1[3] > phase2[3] ? phase1 : phase2
        # Model outputs mass fraction, convert to mg/L
        sim_furfural_mgL = model_EA_rich[2] * rho_org * 1000.0
        push!(results, (Int(run_num), exp_furfural_mgL, sim_furfural_mgL))
    catch e
        @warn "solve_lle failed for run $run_num: $e"
        continue
    end
end

# Export results
output_dir = "TernaryDiagram_2_5"
isdir(output_dir) || mkdir(output_dir)
CSV.write(joinpath(output_dir,"Furfural_concentration_comparison.csv"), results)
println("Furfural comparison CSV saved")

# Box plot
if nrow(results) == 0
    @warn "No results to plot. Exiting boxplot step."
else
    grouped_data = DataFrame(
        Label = repeat(["Experimental", "Simulated"], inner = nrow(results)),
        Furfural_mgL = vcat(results.Exp_FUR_mgL, results.Sim_FUR_mgL)
    )

    grouped_data = dropmissing(grouped_data)

    # make plot
    p = plot(
        ylabel="Furfural in organic phase (mg/L)",
        title="Experimental vs Simulated Furfural in EA-rich Phase",
        legend=:topleft,
        size=(600,400),
        box=:on
    )

    # Box plot
    @df grouped_data boxplot!(:Label, :Furfural_mgL,
        fillalpha=0.3,
        linewidth=1.2,
        outliers=true,
        label=""
    )

    # add jittered points
    for (i, label) in enumerate(["Experimental", "Simulated"])
        data = grouped_data[grouped_data.Label .== label, :Furfural_mgL]
        x = fill(i, length(data)) .+ (rand(length(data)) .- 0.5) .* 0.2  # jitter
        scatter!(x, data,
            markersize=4,
            markerstrokewidth=0,
            alpha=0.6,
            label=""
        )
    end

    # Add mean markers centered over each box
    mean_exp = mean(grouped_data.Furfural_mgL[grouped_data.Label .== "Experimental"])
    mean_sim = mean(grouped_data.Furfural_mgL[grouped_data.Label .== "Simulated"])
    scatter!([1, 2], [mean_exp, mean_sim],
        marker=:diamond,
        color=:red,
        markersize=8,
        label="Mean",
        hover=[@sprintf("%.2f", mean_exp), @sprintf("%.2f", mean_sim)]
    )

    # Add mean val annotations
    annotate!([(1, mean_exp + maximum(grouped_data.Furfural_mgL)*0.02, text(@sprintf("%.0f", mean_exp), 8, :center)),
               (2, mean_sim + maximum(grouped_data.Furfural_mgL)*0.02, text(@sprintf("%.0f", mean_sim), 8, :center))])

    # Save to png
    savefig(joinpath(output_dir,"Furfural_boxplot.png"))
    println("Box plot saved as '$(output_dir)/Furfural_boxplot.png'")
end
