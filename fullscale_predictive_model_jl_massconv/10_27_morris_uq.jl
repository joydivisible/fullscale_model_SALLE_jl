using LinearAlgebra
using CSV, DataFrames, Statistics, Random, Distributions

include("calculate_phase_split.jl")

"""
Run uncertainty propagation with optional Morris sensitivity weighting (micro-UQ)
"""
function uncertainty_analysis(
    tau_file::String,
    feed_compositions::Vector,
    n_samples::Int=20,
    output_dir::String="UncertaintyOutcomes";
    morris_mu_star::Dict = Dict((1,3)=>0.45,(3,1)=>0.20)
)

    println("="^70)
    println("Fast Uncertainty Propagation (Morris-guided micro-UQ)")
    println("="^70)
    println("Tau file: $tau_file")
    println("Feeds: $(length(feed_compositions)), Samples per feed: $n_samples")
    println("="^70)

    if !isdir(output_dir)
        mkdir(output_dir)
        println("Created directory: $output_dir")
    end

    # === Load τ ===
    base_tau = read_tau_matrix(tau_file)

    # Select top 2 most influential τ entries
    sorted = sort(collect(morris_mu_star), by=x->-x[2])
    top_entries = first.(sorted[1:min(2,length(sorted))])
    println("Perturbing τ entries: ", top_entries)
    println()

    all_results = []

    for (feed_idx, z_feed) in enumerate(feed_compositions)
        println("\n", "="^70)
        println("Feed Composition $feed_idx: ", z_feed)
        println("="^70)

        feed_results = []
        for i in 1:n_samples
            print("  Sample $i/$n_samples... ")

            τp = copy(base_tau)
            # perturb only influential τ entries
            for (ii, jj) in top_entries
                τp[ii,jj] *= (1 + 0.02 * (2rand() - 1)) # ±2% random
            end

            # Solve
            try
                x1, x2 = solve_lle(z_feed, τp)

                if any(isnan, x1) || any(isnan, x2)
                    println("NaN")
                    status = "failed"
                    sep = NaN
                else
                    sep = norm(x1 - x2)
                    status = if sep < 0.8 || sep > 1.5 "outlier" else "valid" end
                    println(status == "valid" ? "✓" : "⚠ $(status)")
                end

                push!(feed_results, Dict(
                    :feed_idx => feed_idx,
                    :feed_water => z_feed[1],
                    :feed_furfural => z_feed[2],
                    :feed_ea => z_feed[3],
                    :sample => i,
                    :phase1_water => x1[1],
                    :phase1_furfural => x1[2],
                    :phase1_ea => x1[3],
                    :phase2_water => x2[1],
                    :phase2_furfural => x2[2],
                    :phase2_ea => x2[3],
                    :separation => sep,
                    :status => status
                ))
            catch e
                println("ERROR: $e")
                push!(feed_results, Dict(
                    :feed_idx => feed_idx,
                    :feed_water => z_feed[1],
                    :feed_furfural => z_feed[2],
                    :feed_ea => z_feed[3],
                    :sample => i,
                    :phase1_water => NaN,
                    :phase1_furfural => NaN,
                    :phase1_ea => NaN,
                    :phase2_water => NaN,
                    :phase2_furfural => NaN,
                    :phase2_ea => NaN,
                    :separation => NaN,
                    :status => "error"
                ))
            end
        end

        # === Feed-level stats ===
        feed_df = DataFrame(feed_results)
        valid_df = feed_df[feed_df.status .== "valid", :]
        println("  Valid: $(nrow(valid_df)) / $n_samples")

        if nrow(valid_df) > 0
            println("  Δphase separation = $(round(std(valid_df.separation),digits=3))")
            println("  EA extract = $(round(mean(valid_df.phase2_ea)*100,digits=1)) ± $(round(std(valid_df.phase2_ea)*100,digits=1)) %")
        end

        # Save feed results
        feed_filename = joinpath(output_dir, "feed_$(feed_idx).csv")
        CSV.write(feed_filename, feed_df)
        append!(all_results, feed_results)
    end

    # === Combined stats ===
    println("\n", "="^70)
    println("OVERALL SUMMARY")
    println("="^70)

    all_df = DataFrame(all_results)
    valid_all = all_df[all_df.status .== "valid", :]

    if nrow(valid_all) > 0
        summary = combine(groupby(valid_all, :feed_idx),
            :separation => mean => :sep_mean,
            :separation => std => :sep_std,
            :phase1_water => mean => :water_mean,
            :phase2_ea => mean => :ea_mean,
            :phase2_ea => std => :ea_std)
        CSV.write(joinpath(output_dir, "summary_microUQ.csv"), summary)
        println("Saved summary_microUQ.csv")
    end

    return all_df, valid_all
end


# === Run example ===
feeds = [
    [0.205, 0.010, 0.785],
    [0.308, 0.026, 0.667],
    [0.410, 0.051, 0.538],
    [0.513, 0.103, 0.385],
    [0.226, 0.077, 0.697]
]

all_results, valid_results = uncertainty_analysis(
    "fitted_tau_matrix_2_5.csv",
    feeds,
    20
)
