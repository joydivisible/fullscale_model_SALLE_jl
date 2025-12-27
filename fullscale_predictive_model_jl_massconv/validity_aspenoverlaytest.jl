# Include the necessary helper scripts
include("calculate_phase_split.jl")
include("generate_training_ternary_plots.jl")

"""
    calculate_metrics(exp_points, model_points)

Calculates RMSE and R-squared between experimental and model-predicted compositions.
"""
function main()
    # --- Setup ---
    output_dir = "ValidationDiagram_TrainingSets"
    if !isdir(output_dir)
        mkdir(output_dir)
        println("Created directory: $output_dir")
    end

    tau_files = [
        "fitted_tau_matrix_0.csv", "fitted_tau_matrix_1_75.csv",
        "fitted_tau_matrix_2_5.csv", "fitted_tau_matrix_3_75.csv",
        "fitted_tau_matrix_5.csv", "fitted_tau_matrix_7_5.csv",
        "fitted_tau_matrix_10.csv"
    ]

    data_files = [
        "DATA_SALT_0.csv", "DATA_SALT_1_75.csv",
        "DATA_SALT_2_5.csv", "DATA_SALT_3_75.csv",
        "DATA_SALT_5.csv", "DATA_SALT_7_5.csv",
        "DATA_SALT_10.csv"
    ]

    println("\nStarting batch generation of validation plots...")
    println("-"^50)

    # --- Main Loop ---
    for (i, tau_file) in enumerate(tau_files)
        data_file = data_files[i]
        
        println("Processing: $tau_file with $data_file")

        # --- 1. Load Data and Model ---
        τ_matrix = read_tau_matrix(tau_file)
        aspen_df = CSV.read(joinpath("TrainingSets", data_file), DataFrame)

        # --- 2. Select Tie Lines, Create Intermediate Feeds, and Generate Predictions ---
        total_flows = aspen_df.FLOWWA .+ aspen_df.FLOWEA .+ aspen_df.FLOWFUR
        z_feeds = [[w, f, e] for (w, f, e) in zip(aspen_df.FLOWWA ./ total_flows, aspen_df.FLOWFUR ./ total_flows, aspen_df.FLOWEA ./ total_flows)]
        
        num_tie_lines = length(z_feeds)
        indices_to_plot = round.(Int, range(1, num_tie_lines - 1, length=6)) |> unique # -1 to avoid out of bounds for interpolation

        println("  Selected $(length(indices_to_plot)) base tie lines for display and interpolation.")

        # Store points for different purposes
        model_points_for_metrics = []
        exp_points_for_metrics = []
        model_tie_lines_for_plot = []
        exp_tie_lines_for_plot = []

        exp_p1_water = aspen_df.XWIN1; exp_p1_fur = aspen_df.XFIN1; exp_p1_ea = aspen_df.XEAIN1
        exp_p2_water = aspen_df.XWIN2; exp_p2_fur = aspen_df.XFIN2; exp_p2_ea = aspen_df.XEAIN2

        for idx in indices_to_plot
            # A. Get experimental data for plotting and metrics
            exp_p1 = [exp_p1_water[idx], exp_p1_fur[idx], exp_p1_ea[idx]]
            exp_p2 = [exp_p2_water[idx], exp_p2_fur[idx], exp_p2_ea[idx]]
            push!(exp_tie_lines_for_plot, (exp_p1, exp_p2))
            push!(exp_points_for_metrics, exp_p1, exp_p2)

            # B. Get model prediction for the original feed for fair metrics
            model_p1_metric, model_p2_metric = solve_lle(z_feeds[idx], τ_matrix)
            if !any(isnan, model_p1_metric)
                push!(model_points_for_metrics, model_p1_metric, model_p2_metric)
            end

            # C. Create intermediate feed for plotting model's generalization
            z1 = z_feeds[idx]
            z2 = z_feeds[idx + 1] # The next point in the dataset
            z_intermediate = (z1 .+ z2) ./ 2
            
            # D. Get model prediction for the intermediate feed for plotting
            model_p1_plot, model_p2_plot = solve_lle(z_intermediate, τ_matrix)
            if !any(isnan, model_p1_plot)
                push!(model_tie_lines_for_plot, (model_p1_plot, model_p2_plot))
            end
        end

        # --- 4. Create Plot ---
        plt = plot(size=(800, 750), dpi=300)
        draw_ternary_triangle_clean!()
        add_ternary_grid!(0.1)

        # Separate tie line endpoints by phase type (water-rich vs EA-rich)
        exp_water_rich_cart = []
        exp_ea_rich_cart = []
        model_water_rich_cart = []
        model_ea_rich_cart = []

        # Plot Experimental Tie Lines and identify phases
        for (i_plot, (p1, p2)) in enumerate(exp_tie_lines_for_plot)
            cart_p1 = ternary_to_cartesian(p1[1], p1[3], p1[2]) # W, EA, F
            cart_p2 = ternary_to_cartesian(p2[1], p2[3], p2[2]) # W, EA, F
            
            # Identify which phase is water-rich vs EA-rich
            if p1[1] > p2[1]  # p1 has more water, so p1 is water-rich
                push!(exp_water_rich_cart, cart_p1)
                push!(exp_ea_rich_cart, cart_p2)
            else  # p2 has more water, so p2 is water-rich
                push!(exp_water_rich_cart, cart_p2)
                push!(exp_ea_rich_cart, cart_p1)
            end
            
            # Plot tie line
            plot!([cart_p1[1], cart_p2[1]], [cart_p1[2], cart_p2[2]], 
                  color=:blue, linewidth=2, label=(i_plot==1 ? "Experimental" : ""))
        end

        # Plot Model Tie Lines and identify phases
        for (i_plot, (p1, p2)) in enumerate(model_tie_lines_for_plot)
            cart_p1 = ternary_to_cartesian(p1[1], p1[3], p1[2]) # W, EA, F
            cart_p2 = ternary_to_cartesian(p2[1], p2[3], p2[2]) # W, EA, F
            
            # Identify which phase is water-rich vs EA-rich
            if p1[1] > p2[1]  # p1 has more water, so p1 is water-rich
                push!(model_water_rich_cart, cart_p1)
                push!(model_ea_rich_cart, cart_p2)
            else  # p2 has more water, so p2 is water-rich
                push!(model_water_rich_cart, cart_p2)
                push!(model_ea_rich_cart, cart_p1)
            end
            
            # Plot tie line
            plot!([cart_p1[1], cart_p2[1]], [cart_p1[2], cart_p2[2]], 
                  color=:red, linestyle=:dash, linewidth=2, 
                  label=(i_plot==1 ? "Model (Intermediate Feeds)" : ""))
        end

        # Add markers for experimental tie line endpoints
        if !isempty(exp_water_rich_cart)
            # Water-rich phase - black squares
            scatter!([p[1] for p in exp_water_rich_cart], [p[2] for p in exp_water_rich_cart],
                    markershape=:square, markersize=8, markercolor=:black,
                    markerstrokecolor=:blue, markerstrokewidth=2,
                    label="Exp. Water-rich")
            
            # EA-rich phase - white squares  
            scatter!([p[1] for p in exp_ea_rich_cart], [p[2] for p in exp_ea_rich_cart],
                    markershape=:square, markersize=8, markercolor=:white,
                    markerstrokecolor=:blue, markerstrokewidth=2,
                    label="Exp. EA-rich")
        end

        # Add markers for model tie line endpoints
        if !isempty(model_water_rich_cart)
            # Water-rich phase - black circles
            scatter!([p[1] for p in model_water_rich_cart], [p[2] for p in model_water_rich_cart],
                    markershape=:circle, markersize=8, markercolor=:black,
                    markerstrokecolor=:red, markerstrokewidth=2,
                    label="Model Water-rich")
            
            # EA-rich phase - white circles
            scatter!([p[1] for p in model_ea_rich_cart], [p[2] for p in model_ea_rich_cart],
                    markershape=:circle, markersize=8, markercolor=:white,
                    markerstrokecolor=:red, markerstrokewidth=2,
                    label="Model EA-rich")
        end
        
        plot_title = replace(data_file, "DATA_SALT_" => "", ".csv" => "% NaCl")
        plot!(title="Validation for $plot_title Data", titlefontsize=14,
              legend=:topright, legendfontsize=9)

        # --- 5. Save Plot ---
        plot_name = "validation_plot_" * splitext(data_file)[1] * ".png"
        savefig(plt, joinpath(output_dir, plot_name))
        println("  Plot saved as: $(joinpath(output_dir, plot_name))\n")
    end # End of main for loop

    println("-"^50)
    println("All validation plots generated successfully!")
end # End of main function

# Run the main function
main()