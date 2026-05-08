using CSV, DataFrames, Statistics

# Filepaths
hplc_file = "pre-processed out data.csv"  # HPLC measurements
prep_file = "pre-processed in data.csv"   # Run specifications with Hypo Vol
output_csv = "parsed_Furfural_run.csv"

# Pure-component densities (g/mL)
rho_water, rho_furf, rho_ea = 0.997, 1.159, 0.902
rho_water_gL, rho_furf_gL, rho_ea_gL = rho_water*1000, rho_furf*1000, rho_ea*1000

# Load CSVs
hplc = CSV.read(hplc_file, DataFrame)
prep = CSV.read(prep_file, DataFrame)

# Helper to parse triplicates from HPLC
function parse_triplicates(s)
    mean(parse.(Float64, split(s, '/')))
end

results = []

for i in 1:nrow(prep)
    run = prep[i, Symbol("Run #")]

    # Use Hypo Vol columns directly (mL)
    V_H2O = prep[i, Symbol("Hypo Vol Water (mL)")] / 1000.0  # L
    V_FURF = prep[i, Symbol("Hypo Vol Furfural (mL)")] / 1000.0
    V_EA = prep[i, Symbol("Hypo Vol EA (mL)")] / 1000.0
    V_tot = prep[i, Symbol("Total Vol (mL)")] / 1000.0

    # Masses from exp. vols and densities
    mH2O = V_H2O * rho_water_gL
    mFURF = V_FURF * rho_furf_gL
    mEA = V_EA * rho_ea_gL
    M_tot = mH2O + mFURF + mEA

    # Compute mass fractions
    pct(m) = 100 * m / M_tot

    # HPLC concentrations (mg/L to g/L)
    C_F_a = parse_triplicates(hplc[i, Symbol("Furfural in Water-rich")]) / 1000
    C_EA_a = parse_triplicates(hplc[i, Symbol("Ethyl Acetate in Water-rich")]) / 1000
    C_F_o = parse_triplicates(hplc[i, Symbol("Furfural in EA-rich")]) / 1000
    C_EA_o = parse_triplicates(hplc[i, Symbol("Ethyl Acetate in EA-rich")]) / 1000

    push!(results, (
        run = run,
        V_H2O_L = V_H2O,
        V_FURF_L = V_FURF,
        V_EA_L = V_EA,
        M_H2O_g = mH2O,
        M_FURF_g = mFURF,
        M_EA_g = mEA,
        pct_H2O = pct(mH2O),
        pct_FURF = pct(mFURF),
        pct_EA = pct(mEA),
        C_F_a_gL = C_F_a,
        C_EA_a_gL = C_EA_a,
        C_F_o_gL = C_F_o,
        C_EA_o_gL = C_EA_o
    ))
end

# Export results
df_results = DataFrame(results)
CSV.write(output_csv, df_results)
println("Saved back-calibrated mass% results to $output_csv")
