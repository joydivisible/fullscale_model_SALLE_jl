using CSV, DataFrames, LinearAlgebra, Statistics
include("NRTL_model.jl")  # added to load solve_lle definition
include("calculate_phase_split.jl") 

println("=== TESTING NRTL SOLVER IN YOUR DATA RANGE ===")

# Load your actual data and parameters
println("\n1. LOADING YOUR ACTUAL DATA")
println("="^40)

# NEW: Explicit file checks (instead of try/catch)
script_dir = @__DIR__
tau_candidates = [
    joinpath(script_dir, "fitted_tau_matrix_0.csv"),
    joinpath(pwd(), "fitted_tau_matrix_0.csv"),
    "fitted_tau_matrix_0.csv"
]
train_candidates = [
    joinpath(script_dir, "TrainingSets", "DATA_SALT_0.csv"),
    joinpath(pwd(), "TrainingSets", "DATA_SALT_0.csv"),
    "TrainingSets/DATA_SALT_0.csv"
]

function find_existing(paths)
    for p in paths
        if isfile(p)
            return p
        end
    end
    return nothing
end

tau_path = find_existing(tau_candidates)
train_path = find_existing(train_candidates)

if tau_path === nothing
    println("❌ fitted_tau_matrix_0.csv not found. Searched:")
    for p in tau_candidates
        println("   - ", p)
    end
    println("Current working directory: ", pwd())
    exit(1)
end

if train_path === nothing
    println("❌ TrainingSets/DATA_SALT_0.csv not found. Searched:")
    for p in train_candidates
        println("   - ", p)
    end
    println("Current working directory: ", pwd())
    exit(1)
end

tau_df = CSV.read(tau_path, DataFrame)
# Remove non-numeric column if present
if :Parameter in names(tau_df)
    tau_df = select(tau_df, Not(:Parameter))
end
tau_matrix = Matrix(tau_df)
println("✅ Loaded fitted tau matrix from: ", tau_path, "  size=", size(tau_matrix))
display(round.(tau_matrix, digits=4))

train_df = CSV.read(train_path, DataFrame)
println("✅ Loaded training data from: ", train_path, "  tie-lines=", nrow(train_df))

# Get actual feed compositions from your training data
flows = Matrix(train_df[:, [:FLOWEA, :FLOWWA, :FLOWFUR]])
feed_comps = flows ./ sum(flows, dims=2)

# Get experimental phase compositions
x1_exp_full = Matrix(train_df[:, [:XEAIN1, :XWIN1, :XFIN1]])
x2_exp_full = Matrix(train_df[:, [:XEAIN2, :XWIN2, :XFIN2]])

println("\nData ranges in your training set:")
println("Feed compositions:")
for i in 1:3
    comp_names = ["EA", "WA", "FUR"]
    println("  $(comp_names[i]): $(round(minimum(feed_comps[:, i]), digits=4)) to $(round(maximum(feed_comps[:, i]), digits=4))")
end

println("Experimental phase 1 compositions:")
for i in 1:3
    comp_names = ["EA", "WA", "FUR"]
    println("  $(comp_names[i]): $(round(minimum(x1_exp_full[:, i]), digits=4)) to $(round(maximum(x1_exp_full[:, i]), digits=4))")
end

println("Experimental phase 2 compositions:")
for i in 1:3
    comp_names = ["EA", "WA", "FUR"]
    println("  $(comp_names[i]): $(round(minimum(x2_exp_full[:, i]), digits=4)) to $(round(maximum(x2_exp_full[:, i]), digits=4))")
end

# Test 2: Test solver on actual training data points
println("\n2. TESTING SOLVER ON ACTUAL TRAINING FEEDS")
println("="^50)

# Test on first 5 training points to see if solver can reproduce training data
test_indices = [1, 10, 50, 100, min(200, nrow(train_df))]
test_results = []

for (i, idx) in enumerate(test_indices)
    if idx > nrow(train_df)
        continue
    end
    
    z_test = feed_comps[idx, :]
    x1_exp = x1_exp_full[idx, :]
    x2_exp = x2_exp_full[idx, :]
    
    println("\nTest $i (training point $idx):")
    println("  Feed: z = $(round.(z_test, digits=4))")
    println("  Expected x1 = $(round.(x1_exp, digits=4))")
    println("  Expected x2 = $(round.(x2_exp, digits=4))")
    
    x1_pred, x2_pred = solve_lle(z_test, tau_matrix)
    
    println("  Predicted x1 = $(round.(x1_pred, digits=4))")
    println("  Predicted x2 = $(round.(x2_pred, digits=4))")
    
    # Check sums
    sum1, sum2 = sum(x1_pred), sum(x2_pred)
    println("  Sums: x1=$(round(sum1, digits=4)), x2=$(round(sum2, digits=4))")
    
    # Calculate errors
    error_x1 = norm(x1_pred - x1_exp)
    error_x2 = norm(x2_pred - x2_exp)
    println("  Errors: ||x1_pred - x1_exp|| = $(round(error_x1, digits=4))")
    println("         ||x2_pred - x2_exp|| = $(round(error_x2, digits=4))")
    
    # Store results
    push!(test_results, (
        idx = idx,
        z = z_test,
        x1_exp = x1_exp, x2_exp = x2_exp,
        x1_pred = x1_pred, x2_pred = x2_pred,
        error_x1 = error_x1, error_x2 = error_x2,
        sum1 = sum1, sum2 = sum2
    ))
    
    # Check for major issues
    if abs(sum1 - 1.0) > 0.01 || abs(sum2 - 1.0) > 0.01
        println("  ⚠️ Sum constraint violated!")
    end
    if any(x1_pred .< -0.01) || any(x2_pred .< -0.01)
        println("  ⚠️ Negative mole fractions!")
    end
    if any(x1_pred .> 1.01) || any(x2_pred .> 1.01)
        println("  ⚠️ Mole fractions > 1!")
    end
    if error_x1 > 0.1 || error_x2 > 0.1
        println("  ⚠️ Large prediction error!")
    end
end

# Test 3: Summary of results
println("\n3. SUMMARY OF SOLVER PERFORMANCE")
println("="^45)

successful_tests = [r for r in test_results if r.x1_pred !== nothing]
failed_tests = [r for r in test_results if r.x1_pred === nothing]

println("Successful predictions: $(length(successful_tests))/$(length(test_results))")
println("Failed predictions: $(length(failed_tests))")

if length(successful_tests) > 0
    errors_x1 = [r.error_x1 for r in successful_tests]
    errors_x2 = [r.error_x2 for r in successful_tests]
    
    println("\nError statistics:")
    println("  X1 errors - Mean: $(round(mean(errors_x1), digits=4)), Max: $(round(maximum(errors_x1), digits=4))")
    println("  X2 errors - Mean: $(round(mean(errors_x2), digits=4)), Max: $(round(maximum(errors_x2), digits=4))")
    
    # Check sum constraints
    sum_errors_x1 = [abs(r.sum1 - 1.0) for r in successful_tests]
    sum_errors_x2 = [abs(r.sum2 - 1.0) for r in successful_tests]
    
    println("  Sum constraint violations:")
    println("    X1 - Mean: $(round(mean(sum_errors_x1), digits=6)), Max: $(round(maximum(sum_errors_x1), digits=6))")
    println("    X2 - Mean: $(round(mean(sum_errors_x2), digits=6)), Max: $(round(maximum(sum_errors_x2), digits=6))")
end

# Test 4: Test validation feeds (the ones causing negative R²)
println("\n4. TESTING VALIDATION FEEDS (FROM NEGATIVE R² CASE)")
println("="^55)

# Recreate the exact same validation feeds that caused the issue
Random.seed!(42)
validation_feeds = []
N = size(feed_comps, 1)

for i in 1:6
    i1, i2 = rand(1:N), rand(1:N)
    α = rand()
    z = α * feed_comps[i1, :] + (1 - α) * feed_comps[i2, :]
    push!(validation_feeds, z ./ sum(z))
end

println("Testing the exact validation feeds that caused negative R²:")

validation_results = []
for (i, z) in enumerate(validation_feeds)
    println("\nValidation feed $i: z = $(round.(z, digits=4))")
    
    # Get interpolated experimental values (same as in your validation)
    dists = [norm(z - feed_comps[j, :]) for j in 1:size(feed_comps,1)]
    idxs = sortperm(dists)[1:2]
    i1, i2 = idxs
    d1, d2 = dists[i1], dists[i2]
    
    if d1 + d2 ≈ 0
        w1, w2 = 0.5, 0.5
    else
        w1 = d2 / (d1 + d2)
        w2 = d1 / (d1 + d2)
    end
    
    x1_exp_interp = w1 * x1_exp_full[i1, :] + w2 * x1_exp_full[i2, :]
    x2_exp_interp = w1 * x2_exp_full[i1, :] + w2 * x2_exp_full[i2, :]
    
    println("  Interpolated experimental:")
    println("    x1_exp = $(round.(x1_exp_interp, digits=4))")
    println("    x2_exp = $(round.(x2_exp_interp, digits=4))")
    
    x1_pred, x2_pred = solve_lle(z, tau_matrix)
    println("  Predicted:")
    println("    x1_pred = $(round.(x1_pred, digits=4))")
    println("    x2_pred = $(round.(x2_pred, digits=4))")
    
    # Calculate component-wise errors that go into R²
    errors_x1_comp = (x1_pred - x1_exp_interp).^2
    errors_x2_comp = (x2_pred - x2_exp_interp).^2
    
    println("  Squared errors by component:")
    println("    X1: EA=$(round(errors_x1_comp[1], digits=6)), WA=$(round(errors_x1_comp[2], digits=6)), FUR=$(round(errors_x1_comp[3], digits=6))")
    println("    X2: EA=$(round(errors_x2_comp[1], digits=6)), WA=$(round(errors_x2_comp[2], digits=6)), FUR=$(round(errors_x2_comp[3], digits=6))")
    
    push!(validation_results, (
        z = z, x1_exp = x1_exp_interp, x2_exp = x2_exp_interp,
        x1_pred = x1_pred, x2_pred = x2_pred
    ))
end

# Test 5: Manual R² calculation to understand the negative values
if length(validation_results) > 0
    println("\n5. MANUAL R² CALCULATION")
    println("="^35)
    
    # Collect all predictions and experimental values
    all_x1_pred = hcat([r.x1_pred for r in validation_results]...)'
    all_x2_pred = hcat([r.x2_pred for r in validation_results]...)'
    all_x1_exp = hcat([r.x1_exp for r in validation_results]...)'
    all_x2_exp = hcat([r.x2_exp for r in validation_results]...)'
    
    comp_names = ["EA", "WA", "FUR"]
    
    println("Manual R² calculation for each component:")
    for j in 1:3
        # Phase 1
        y_true_x1 = all_x1_exp[:, j]
        y_pred_x1 = all_x1_pred[:, j]
        
        ss_res_x1 = sum((y_true_x1 - y_pred_x1).^2)
        ss_tot_x1 = sum((y_true_x1 .- mean(y_true_x1)).^2)
        r2_x1 = 1 - ss_res_x1/ss_tot_x1
        
        # Phase 2
        y_true_x2 = all_x2_exp[:, j]
        y_pred_x2 = all_x2_pred[:, j]
        
        ss_res_x2 = sum((y_true_x2 - y_pred_x2).^2)
        ss_tot_x2 = sum((y_true_x2 .- mean(y_true_x2)).^2)
        r2_x2 = 1 - ss_res_x2/ss_tot_x2
        
        println("\n$(comp_names[j]) component:")
        println("  X1: SS_res=$(round(ss_res_x1, digits=6)), SS_tot=$(round(ss_tot_x1, digits=6)), R²=$(round(r2_x1, digits=4))")
        println("  X2: SS_res=$(round(ss_res_x2, digits=6)), SS_tot=$(round(ss_tot_x2, digits=6)), R²=$(round(r2_x2, digits=4))")
        
        println("  X1 range: exp [$(round(minimum(y_true_x1), digits=4)), $(round(maximum(y_true_x1), digits=4))], pred [$(round(minimum(y_pred_x1), digits=4)), $(round(maximum(y_pred_x1), digits=4))]")
        println("  X2 range: exp [$(round(minimum(y_true_x2), digits=4)), $(round(maximum(y_true_x2), digits=4))], pred [$(round(minimum(y_pred_x2), digits=4)), $(round(maximum(y_pred_x2), digits=4))]")
        
        if ss_tot_x1 < 1e-10
            println("  ⚠️ X1: Very small experimental variance - R² calculation unreliable")
        end
        if ss_tot_x2 < 1e-10
            println("  ⚠️ X2: Very small experimental variance - R² calculation unreliable")
        end
    end
end

println("\n6. DIAGNOSTIC CONCLUSIONS")
println("="^35)
println("Run this test to identify:")
println("- Whether your solver works on training data")
println("- What specific errors occur in validation")
println("- Why R² values are negative")
println("- Whether the issue is solver convergence or parameter quality")