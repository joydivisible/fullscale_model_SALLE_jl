# NRTL_Model.jl
# This file contains the definitive, correct, and simplified implementation of the
# NRTL activity coefficient model. It uses explicit for-loops for clarity and robustness.

using LinearAlgebra

# Define the non-randomness constant globally for this model
const ALPHA_CONST = 0.2

"""
    activity_coefficients(x, tau)

Calculates the activity coefficients (γ) for a mixture with mole fractions `x`
using the NRTL model with binary interaction parameters `tau`.

# Arguments
- `x`: A vector of mole fractions for each component.
- `tau`: A matrix of binary interaction parameters (τ_ij).

# Returns
- A vector of activity coefficients (γ) for each component.
"""
function activity_coefficients(x, tau)
    n = length(x)
    G = exp.(-ALPHA_CONST .* tau)
    ε = 1e-8

    S = [sum(x[j] * G[j,i] for j in 1:n) for i in 1:n]

    ln_gamma = [
        begin
            term1 = sum(x[j] * G[j,i] * tau[j,i] for j in 1:n) / (S[i] + ε)
            term2 = sum(
                (x[j] * G[j,i] / (S[j] + ε)) *
                (tau[j,i] - sum(x[k] * G[k,i] * tau[k,i] for k in 1:n) / (S[i] + ε))
                for j in 1:n
            )
            term1 + term2
        end
        for i in 1:n
    ]

    return exp.(ln_gamma)
end
