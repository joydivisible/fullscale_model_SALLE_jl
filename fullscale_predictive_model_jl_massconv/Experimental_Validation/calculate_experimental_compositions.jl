using CSV, DataFrames, Statistics, Interpolations, LinearAlgebra

# INPUT DATA
# Component densities (g/mL)
ρ_water = 1.00
ρ_furfural = 1.16
ρ_EA = 0.902

# Molar masses (g/mol)
M_water = 18.015
M_EA    = 88.106

# SHRINKAGE CORRECTION FROM EXCESS MOLAR VOLUME (not currently used in final calculations, but kept for potential future use)
# Molar volumes (cm³/mol)
V_H2O = M_water / ρ_water
V_EA  = M_EA / ρ_EA

# Approx. VE(xEA) data for water + EA binary (318 K)
VE_data = [
    (0.0,  0.0),
    (0.1, -1.2),
    (0.2, -1.1),
    (0.4, -0.8),
    (0.6, -0.6),
    (0.8, -0.4),
    (1.0,  0.0)
]

xvals = [d[1] for d in VE_data]
yvals = [d[2] for d in VE_data]
VE_interp = LinearInterpolation(xvals, yvals; extrapolation_bc=Line())

function q_shrink(xEA::Float64)
    VE = VE_interp(clamp(xEA, 0, 1))
    Videal = xEA * V_EA + (1 - xEA) * V_H2O
    return 1 + VE / Videal      # <1 means contraction
end

# DATA LOADING
in_df  = CSV.read("pre-processed in data.csv", DataFrame)
out_df = CSV.read("pre-processed out data.csv", DataFrame)

# FEED-BASED SHRINKAGE FACTOR
n_water_feed = in_df."Amount Water (g)" ./ M_water
n_ea_feed    = in_df."Amount Ethyl Acetate (g)" ./ M_EA
xEA_feed     = n_ea_feed ./ (n_ea_feed .+ n_water_feed)
q_factors    = q_shrink.(xEA_feed)

# PHASE VOLUMES
# assumption is that water phase does not shrink/shrinks negligibly -- EA phase does
vol_w = Float64.(in_df."Hypo Vol Water (mL)")
vol_e = Float64.(in_df."Hypo Vol EA (mL)") .* q_factors   # shrinkage-corrected EA volume

# HELPER: triplicate parser
function parse_triplicate_mean(col)
    [ismissing(v) || isempty(String(v)) ? missing :
        mean(parse.(Float64, split(String(v), "/"))) for v in col]
end

# FEED MASSES (g)
m_water_feed = in_df."Amount Water (g)"
m_furf_feed  = in_df."Amount Furfural (g)"
m_ea_feed    = in_df."Amount Ethyl Acetate (g)"

# HPLC DATA (mg/L)
furfural_w = parse_triplicate_mean(out_df."Furfural in Water-rich")
ea_w       = parse_triplicate_mean(out_df."Ethyl Acetate in Water-rich")
furfural_w = coalesce.(furfural_w, 0.0)
ea_w       = coalesce.(ea_w, 0.0)

# CALCULATIONS
# Solute masses in water-rich phase
furf_mass_w = furfural_w .* vol_w ./ 1e6
ea_mass_w   = ea_w .* vol_w ./ 1e6

# Estimate water-rich phase total mass
mass_w_total = vol_w .* ρ_water  # g

# Mass of water in water-rich phase
m_wat_w = mass_w_total .- furf_mass_w .- ea_mass_w

# Normalizing water-rich compositions
total_w   = m_wat_w .+ furf_mass_w .+ ea_mass_w
xw_water  = m_wat_w   ./ total_w
xw_furf   = furf_mass_w ./ total_w
xw_ea     = ea_mass_w   ./ total_w

# get EA-rich phase by mass balance
m_wat_e = m_water_feed .- m_wat_w
m_fur_e = m_furf_feed  .- furf_mass_w
m_ea_e  = m_ea_feed    .- ea_mass_w

# Normalize EA-rich compositions
total_e  = m_wat_e .+ m_fur_e .+ m_ea_e
xe_water = m_wat_e ./ total_e
xe_furf  = m_fur_e ./ total_e
xe_ea    = m_ea_e  ./ total_e

# CREATE OUTPUT DATAFRAME
results_df = DataFrame(
    Run = 1:length(xw_water),
    Water_rich_W  = xw_water,
    Water_rich_F  = xw_furf,
    Water_rich_EA = xw_ea,
    EA_rich_W     = xe_water,
    EA_rich_F     = xe_furf,
    EA_rich_EA    = xe_ea,
    Shrink_factor = q_factors,
    Water_phase_vol_mL = vol_w,
    EA_phase_vol_corrected_mL = vol_e,
    Water_phase_mass_g = mass_w_total,
    EA_phase_mass_g    = total_e,
    Furfural_in_water_g = furf_mass_w,
    EA_in_water_g       = ea_mass_w
)

# PRINT RESULTS
println("WATER-RICH PHASE COMPOSITIONS (W, F, EA)")
for i in 1:nrow(results_df)
    println("Run $i: ($(round(results_df.Water_rich_W[i], digits=4)), ",
        "$(round(results_df.Water_rich_F[i], digits=4)), ",
        "$(round(results_df.Water_rich_EA[i], digits=4)))")
end
println()

println("EA-RICH PHASE COMPOSITIONS (W, F, EA)")
for i in 1:nrow(results_df)
    println("Run $i: ($(round(results_df.EA_rich_W[i], digits=4)), ",
        "$(round(results_df.EA_rich_F[i], digits=4)), ",
        "$(round(results_df.EA_rich_EA[i], digits=4)))")
end
println()

# MASS BALANCE CHECK
println("MASS BALANCE CLOSURE CHECK")
total_feed      = m_water_feed .+ m_furf_feed .+ m_ea_feed
total_recovered = mass_w_total .+ total_e
mass_balance_error = abs.((total_recovered .- total_feed) ./ total_feed) .* 100

println("Run | Feed (g) | Recovered (g) | Error (%)")
for i in 1:length(total_feed)
    println("$i | $(round(total_feed[i], digits=4)) | ",
        "$(round(total_recovered[i], digits=4)) | ",
        "$(round(mass_balance_error[i], digits=2))")
end
println()

# WRITE INTO CSV
CSV.write("ternary_compositions.csv", results_df)
println("Results written to: ternary_compositions.csv")



