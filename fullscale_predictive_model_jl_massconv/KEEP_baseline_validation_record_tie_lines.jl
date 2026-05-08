using CSV
using DataFrames

include("calculate_phase_split.jl")

const MW = [18.01528, 96.084, 88.106]  # H2O, Furfural, EA

function mole_to_mass(x)
    w = x .* MW
    return w ./ sum(w)
end

function mass_to_mole(w)
    x = w ./ MW
    return x ./ sum(x)
end

function classify_phases(p1::Vector{Float64}, p2::Vector{Float64})
    # higher water fraction -> water-rich
    if p1[1] > p2[1]
        return p1, p2   # water-rich, ea-rich
    else
        return p2, p1   # water-rich, ea-rich
    end
end

function main()
    output_dir = "BaselineValidationData"
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

    println("\nGenerating validation datasets with interpolated reference tie lines...")
    println("-"^70)

    for (i, tau_file) in enumerate(tau_files)
        data_file = data_files[i]
        salt_tag = replace(data_file, "DATA_SALT_" => "", ".csv" => "")

        println("Processing: $tau_file with $data_file")

        τ_matrix = read_tau_matrix(tau_file)
        aspen_df = CSV.read(joinpath("TrainingSets", data_file), DataFrame)

        # Feed prep: flow-based mass fractions -> mole fractions
        total_flows = aspen_df.FLOWWA .+ aspen_df.FLOWEA .+ aspen_df.FLOWFUR
        z_feeds_mass = [[w, f, e] for (w, f, e) in zip(
            aspen_df.FLOWWA ./ total_flows,
            aspen_df.FLOWFUR ./ total_flows,
            aspen_df.FLOWEA ./ total_flows
        )]
        z_feeds = [mass_to_mole(z) for z in z_feeds_mass]

        num_tie_lines = length(z_feeds)
        indices_to_plot = round.(Int, range(1, num_tie_lines - 1, length = 6)) |> unique

        exp_df = DataFrame(
            idx = Int[],
            water_rich_w = Float64[],
            water_rich_f = Float64[],
            water_rich_e = Float64[],
            ea_rich_w = Float64[],
            ea_rich_f = Float64[],
            ea_rich_e = Float64[]
        )

        model_df = DataFrame(
            idx = Int[],
            water_rich_w = Float64[],
            water_rich_f = Float64[],
            water_rich_e = Float64[],
            ea_rich_w = Float64[],
            ea_rich_f = Float64[],
            ea_rich_e = Float64[]
        )

        exp_p1_water = aspen_df.XWIN1
        exp_p1_fur   = aspen_df.XFIN1
        exp_p1_ea    = aspen_df.XEAIN1

        exp_p2_water = aspen_df.XWIN2
        exp_p2_fur   = aspen_df.XFIN2
        exp_p2_ea    = aspen_df.XEAIN2

        for idx in indices_to_plot
            # average experimental phase compositions at idx and idx+1, then convert mass -> mole
            exp_p1_mass = Float64[
                (exp_p1_water[idx] + exp_p1_water[idx + 1]) / 2,
                (exp_p1_fur[idx]   + exp_p1_fur[idx + 1])   / 2,
                (exp_p1_ea[idx]    + exp_p1_ea[idx + 1])    / 2
            ]
            exp_p2_mass = Float64[
                (exp_p2_water[idx] + exp_p2_water[idx + 1]) / 2,
                (exp_p2_fur[idx]   + exp_p2_fur[idx + 1])   / 2,
                (exp_p2_ea[idx]    + exp_p2_ea[idx + 1])    / 2
            ]

            exp_p1 = mass_to_mole(exp_p1_mass)
            exp_p2 = mass_to_mole(exp_p2_mass)

            exp_wr, exp_ea = classify_phases(exp_p1, exp_p2)

            push!(exp_df, (
                idx,
                exp_wr[1], exp_wr[2], exp_wr[3],
                exp_ea[1], exp_ea[2], exp_ea[3]
            ))

            # interpolate feed in mole space, then solve
            z1 = z_feeds[idx]
            z2 = z_feeds[idx + 1]
            z_intermediate = (z1 .+ z2) ./ 2

            model_p1, model_p2 = solve_lle(z_intermediate, τ_matrix)

            if !any(isnan, model_p1) && !any(isnan, model_p2)
                model_wr, model_ea = classify_phases(
                    Float64[model_p1[1], model_p1[2], model_p1[3]],
                    Float64[model_p2[1], model_p2[2], model_p2[3]]
                )

                push!(model_df, (
                    idx,
                    model_wr[1], model_wr[2], model_wr[3],
                    model_ea[1], model_ea[2], model_ea[3]
                ))
            end
        end

        exp_file = joinpath(output_dir, "exp_validation_data_$salt_tag.csv")
        model_file = joinpath(output_dir, "model_validation_data_$salt_tag.csv")

        CSV.write(exp_file, exp_df)
        CSV.write(model_file, model_df)

        println("  Saved: $exp_file")
        println("  Saved: $model_file")
        println()
    end

    println("-"^70)
    println("All validation data exported.")
end

main()