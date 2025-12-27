using CSV
using DataFrames
using Plots

# --- CONFIGURATION ---
input_folder = "."
output_folder = "COMPARE_TAU_PLOTS"

# Create output folder if it doesn't exist
isdir(output_folder) || mkpath(output_folder)

# Helper: extract numeric value from filename (e.g., "1_75" → 1.75)
function extract_numeric(fname::String)
    base = splitext(fname)[1]            # remove .csv
    num_str = replace(base, "_" => ".")  # underscores → dots
    m = match(r"([0-9]+\.[0-9]+|[0-9]+)", num_str)
    return m === nothing ? NaN : parse(Float64, m.captures[1])
end

# Get list of CSV files sorted by numeric value in filename
csv_files = sort(
    filter(f -> endswith(lowercase(f), ".csv") &&
                 occursin("tau", lowercase(f)),
           readdir(input_folder)),
    by = extract_numeric
)

# Read matrices into Dict
matrices = Dict{String, Matrix{Float64}}()
numeric_labels = Float64[]  # store numeric values for x-axis

for file in csv_files
    path = joinpath(input_folder, file)
    df = CSV.read(path, DataFrame)           # header row: x1, x2, x3
    df_matrix = Matrix{Float64}(df)          # convert to numeric
    matrices[file] = df_matrix
    push!(numeric_labels, extract_numeric(file))
end

# Off-diagonal indices
off_diag_indices = [(i, j) for i in 1:3 for j in 1:3 if i != j]

# Plot each element
for (i, j) in off_diag_indices
    y_values = [matrices[fname][i, j] for fname in csv_files]

    plt = plot(
        numeric_labels, y_values, marker=:o,
        xlabel="Salt amount [% (m/m)]", ylabel="Value",
        title="Tau matrix element ($i,$j)",
        legend=false, grid=true
    )

    savefig(plt, joinpath(output_folder, "COMPARE_TAU_$(i),$(j).png"))
end

println("Plots saved to '$output_folder'")

