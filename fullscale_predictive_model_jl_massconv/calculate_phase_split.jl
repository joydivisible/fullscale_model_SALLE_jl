#=
This script calculates the LLE phase split.

=#

using LinearAlgebra
using Optimization, OptimizationOptimJL
using CSV
using DataFrames
using Statistics
import Zygote

include("NRTL_Model.jl")

# --- Helper Functions ---

function read_tau_matrix(filepath::String)
    df = CSV.read(filepath, DataFrame)
    tau_matrix = Matrix(df)
    println("Loaded tau matrix from $(filepath):")
    display(tau_matrix)
    return tau_matrix
end

"""
    lle_objective(x_guess, z, tau)

A consistent and AD-friendly objective function for LLE.
The error term `(ln(γx₁)-ln(γx₂))²` matches the training script's loss function.
"""
function lle_objective(x_guess, z, tau)
    n = length(z)
    # Both phases are optimization variables
    x = @view x_guess[1:n]
    y = @view x_guess[n+1:end]
    
    ε = 1e-12 # Use a small epsilon for numerical stability
    
    x_sum = sum(x)
    y_sum = sum(y)

    # Normalize compositions safely
    x_norm = x ./ (x_sum + ε)
    y_norm = y ./ (y_sum + ε)

    gamma_x = activity_coefficients(x_norm, tau)
    gamma_y = activity_coefficients(y_norm, tau)

    # 1. Phase Equilibrium Error (consistent with training script)
    # Using log difference of activities, adding epsilon to prevent log(0)
    log_ax = log.(x_norm .* gamma_x .+ ε)
    log_ay = log.(y_norm .* gamma_y .+ ε)
    phase_eq_error = sum((log_ax .- log_ay).^2)

    # 2. Mass Balance Error
    diff_xy = x_norm .- y_norm
    denom = sum(diff_xy .* diff_xy) + ε
    β_opt = sum((z .- y_norm) .* diff_xy) / denom
    z_calc = β_opt .* x_norm .+ (1 - β_opt) .* y_norm
    mass_balance_error = 10.0 * sum((z .- z_calc).^2) # Weighted more heavily

    # 3. Differentiable Penalty for β being outside [0, 1]
    beta_penalty = 100.0 * (max(0.0, -β_opt)^2 + max(0.0, β_opt - 1.0)^2)

    # 4. Differentiable Penalty for trivial solutions
    phase_separation = sum((x_norm .- y_norm).^2)
    trivial_penalty = 1000.0 * exp(-50.0 * phase_separation) # Made steeper

    # Combine all errors/penalties
    total_error = phase_eq_error + mass_balance_error + beta_penalty + trivial_penalty

    return isfinite(total_error) ? total_error : 1e6
end


function solve_lle(z, tau)
    n = length(z)
    z_norm = z ./ sum(z)
    best_x_global, best_y_global = nothing, nothing
    best_loss_global = Inf

    # More diverse initial strategies
    initial_strategies = [
        (z_norm .+ 0.2, z_norm .- 0.2), # Push apart
        (z_norm .* 0.1, z_norm .* 1.9), # Stretch
        ([0.8, 0.1, 0.1], [0.1, 0.1, 0.8]), # Water-rich vs EA-rich
        ([0.1, 0.1, 0.8], [0.8, 0.1, 0.1]), # EA-rich vs Water-rich
        ([0.5, 0.5, 0.01], [0.01, 0.01, 0.98]) # High separation guesses
    ]

    for (i, (x_init, y_init)) in enumerate(initial_strategies)
        x0 = vcat(x_init, y_init)
        x0 = clamp.(x0, 1e-6, 1.0 - 1e-6)
        
        obj_func = OptimizationFunction((x, p) -> lle_objective(x, z_norm, tau), Optimization.AutoZygote())
        prob = OptimizationProblem(obj_func, x0, nothing, lb = fill(1e-6, length(x0)), ub = fill(1.0, length(x0)))
        
        # Use a robust optimizer with more iterations
        sol = solve(prob, OptimizationOptimJL.LBFGS(), maxiters = 300, reltol=1e-8)
        
        if sol.retcode == ReturnCode.Success || sol.retcode == ReturnCode.MaxIters
            final_loss = lle_objective(sol.u, z_norm, tau)
            if final_loss < best_loss_global
                best_loss_global = final_loss
                best_x_global = sol.u[1:n]
                best_y_global = sol.u[n+1:end]
            end
        end
    end

    if best_x_global !== nothing
        x_eq = best_x_global ./ sum(best_x_global)
        y_eq = best_y_global ./ sum(best_y_global)
    else
        x_eq, y_eq = (fill(NaN, n), fill(NaN, n))
        println("Warning: Solver failed to find a solution.")
    end

    # Return phases sorted by water content for consistency
    if x_eq[1] < y_eq[1]
        return y_eq, x_eq
    else
        return x_eq, y_eq
    end
end


# --- Main Execution Block ---

function main()
    τ_fitted = read_tau_matrix("fitted_tau_matrix_0.csv")

    # CRITICAL: This vector MUST follow the component order used for fitting:
    # Index 1: Water, Index 2: Furfural, Index 3: Ethyl Acetate
    z_feed = [0.2, 0.01, 0.79]

    println("\nCalculating LLE split with the following feed composition:")
    println("Feed Composition (z): [Water, Furfural, EA]")
    display(z_feed)
    println("-"^40)

    x1_eq, x2_eq = solve_lle(z_feed, τ_fitted)

    println("\nCalculation Complete.")
    println("Predicted Equilibrium Phases (Water-rich, EA-rich):")
    println("  Phase 1: ", round.(x1_eq, digits=4))
    println("  Phase 2: ", round.(x2_eq, digits=4))
    println("-"^40)
end

