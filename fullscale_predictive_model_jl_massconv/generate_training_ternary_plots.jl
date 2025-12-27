using CSV, DataFrames, Plots

println("Loading data and creating ternary diagram...")


# Function to convert ternary coordinates to cartesian for plotting
function ternary_to_cartesian(a, b, c)
    # Normalize to ensure a + b + c = 1
    total = a + b + c
    a, b, c = a/total, b/total, c/total
    
    # Convert to cartesian coordinates
    x = 0.5 * (2*b + c)
    y = (sqrt(3)/2) * c
    return x, y
end


using CSV, DataFrames, Plots

println("Loading data and creating ternary diagram...")

# Create output directory if it doesn't exist
output_dir = "TernaryDiagram_TrainingSets"
if !isdir(output_dir)
    mkdir(output_dir)
end

# Function to convert ternary coordinates to cartesian for plotting
function ternary_to_cartesian(a, b, c)
    # Normalize to ensure a + b + c = 1
    total = a + b + c
    a, b, c = a/total, b/total, c/total
    
    # Convert to cartesian coordinates
    x = 0.5 * (2*b + c)
    y = (sqrt(3)/2) * c
    return x, y
end

# Function to add oriented axis ticks
function add_axis_ticks!()
    tick_values = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    tick_length = 0.015
    
    for val in tick_values
        # Bottom axis ticks (Water to EA axis) - horizontal ticks
        x_pos = val
        y_pos = 0.0
        plot!([x_pos, x_pos], [y_pos - tick_length, y_pos + tick_length], 
              color=:black, linewidth=1, label="")
        
        # Left axis ticks (Water to Furfural) - angled at 60°
        # Points along left edge: from (0,0) to (0.5, sqrt(3)/2)
        x_left = 0.5 * val
        y_left = (sqrt(3)/2) * val
        # Perpendicular to left edge (60° rotation from horizontal)
        tick_dx = tick_length * cos(π/6)  # cos(30°) for perpendicular to 60° line
        tick_dy = tick_length * sin(π/6)  # sin(30°)
        plot!([x_left - tick_dx, x_left + tick_dx], 
              [y_left + tick_dy, y_left - tick_dy], 
              color=:black, linewidth=1, label="")
        
        # Right axis ticks (EA to Furfural) - angled at -60°
        # Points along right edge: from (1,0) to (0.5, sqrt(3)/2)  
        x_right = 1.0 - 0.5 * val
        y_right = (sqrt(3)/2) * val
        # Perpendicular to right edge (-60° rotation from horizontal)
        tick_dx = tick_length * cos(π/6)  # cos(30°) 
        tick_dy = tick_length * sin(π/6)  # sin(30°)
        plot!([x_right + tick_dx, x_right - tick_dx], 
              [y_right + tick_dy, y_right - tick_dy], 
              color=:black, linewidth=1, label="")
    end
    
    # Add tick labels with proper orientation
    add_tick_labels!()
end


# Fixed tick labels function
function add_tick_labels!()
    tick_values = [0.2, 0.4, 0.6, 0.8]

    # Left edge: Water to Furfural (uphill, shift closer to Furfural by half a tick)
    tick_shift = 0.05  # half a tick for tick_values = [0.2, 0.4, 0.6, 0.8]
    for (i, val) in enumerate(tick_values)
        shifted_val = val + tick_shift
        x_left = 0.5 * shifted_val - 0.06
        y_left = (sqrt(3)/2) * shifted_val
        annotate!(x_left, y_left, text(string(round(val, sigdigits=2)), 10, :black, rotation=60, halign=:right, valign=:center))
    end

    # Right edge: Furfural to EA (downhill, shift closer to Furfural by half a tick)
    for (i, val) in enumerate(tick_values)
        rev_val = 1.0 - val
        shifted_rev_val = rev_val + tick_shift
        x_right = 1.0 - 0.5 * shifted_rev_val + 0.06
        y_right = (sqrt(3)/2) * shifted_rev_val
        annotate!(x_right, y_right, text(string(round(val, sigdigits=2)), 10, :black, rotation=-60, halign=:left, valign=:center))
    end

    # Base: EA to Water (leftward)
    for val in tick_values
        x_base = 1.0 - val
        annotate!(x_base, -0.05, text(string(round(val, sigdigits=2)), 10, :black, rotation=0, halign=:center, valign=:top))
    end

    # Vertices labeled as "1.0"
    annotate!(-0.06, -0.02, text("1.0", 10, :black, rotation=0, halign=:right, valign=:center))
    annotate!(1.06, -0.02, text("1.0", 10, :black, rotation=0, halign=:left, valign=:center))
    annotate!(0.5, sqrt(3)/2 + 0.04, text("1.0", 10, :black, rotation=0, halign=:center, valign=:bottom))
end

# Function to draw ternary triangle with properly oriented axis labels and ticks
function draw_ternary_triangle_clean!()
    # Triangle vertices
    vertices_x = [0.0, 1.0, 0.5, 0.0]
    vertices_y = [0.0, 0.0, sqrt(3)/2, 0.0]
    plot!(vertices_x, vertices_y,
          color=:black, linewidth=2, label="",
          aspect_ratio=:equal,
          showaxis=false,
          grid=false,
          xlims=(-0.15, 1.15),
          ylims=(-0.15, 1.0))

    # Axis labels positioned slightly away from vertices, aligned with each axis
    # Water: bottom-left, label positioned below and slightly left
    annotate!(-0.05, -0.08, text("Water", 14, :black, rotation=0, halign=:center, valign=:top))
    
    # Ethyl Acetate: bottom-right, label positioned below and slightly right  
    annotate!(1.05, -0.08, text("Ethyl Acetate", 14, :black, rotation=0, halign=:center, valign=:top))
    
    # Furfural: top vertex, label positioned above
    annotate!(0.5, sqrt(3)/2 + 0.08, text("Furfural", 14, :black, rotation=0, halign=:center, valign=:bottom))
    
    # Add axis ticks with proper orientation
    add_axis_ticks!()
end

# Function to add grid lines
function add_ternary_grid!(step=0.1)
    for i in step:step:(1-step)
        # Horizontal lines (constant c - Furfural)
        x1, y1 = ternary_to_cartesian(1-i, 0, i)
        x2, y2 = ternary_to_cartesian(0, 1-i, i)
        plot!([x1, x2], [y1, y2], color=:gray, alpha=0.3, linewidth=0.5, label="")
        
        # Left diagonal lines (constant a - Water)
        x1, y1 = ternary_to_cartesian(i, 1-i, 0)
        x2, y2 = ternary_to_cartesian(i, 0, 1-i)
        plot!([x1, x2], [y1, y2], color=:gray, alpha=0.3, linewidth=0.5, label="")
        
        # Right diagonal lines (constant b - EA)
        x1, y1 = ternary_to_cartesian(1-i, i, 0)
        x2, y2 = ternary_to_cartesian(0, i, 1-i)
        plot!([x1, x2], [y1, y2], color=:gray, alpha=0.3, linewidth=0.5, label="")
    end
end

# Get all CSV files in TrainingSets
csv_files = filter(x -> endswith(x, ".csv"), readdir("TrainingSets"))

# Loop through each CSV file
for csv_file in csv_files
    data_file = joinpath("TrainingSets", csv_file)
    println("\\nProcessing $data_file ...")
    
    if !isfile(data_file)
        println("Error: $data_file not found in TrainingSets directory")
        continue
    end
    
    try
        # Load data
        df = CSV.read(data_file, DataFrame)
        println("Loaded $(nrow(df)) rows from $data_file")
        
        # Filter out invalid rows if you have a Status column
        if "Status" in names(df)
            df = df[df.Status .== "OK", :]
            println("After filtering: $(nrow(df)) valid rows")
        end
        
        # Extract compositions using exact column names
        xeain1 = df[:, "XEAIN1"]
        xwin1 = df[:, "XWIN1"]
        xfin1 = df[:, "XFIN1"]
        xeain2 = df[:, "XEAIN2"]
        xwin2 = df[:, "XWIN2"]
        xfin2 = df[:, "XFIN2"]
        
        println("Successfully extracted composition data from specified columns")
        
        # Calculate normalized compositions (ensure they sum to 1)
        EA1 = xeain1 ./ (xeain1 .+ xwin1 .+ xfin1)
        WA1 = xwin1 ./ (xeain1 .+ xwin1 .+ xfin1)
        FUR1 = xfin1 ./ (xeain1 .+ xwin1 .+ xfin1)
        
        EA2 = xeain2 ./ (xeain2 .+ xwin2 .+ xfin2)
        WA2 = xwin2 ./ (xeain2 .+ xwin2 .+ xfin2)
        FUR2 = xfin2 ./ (xeain2 .+ xwin2 .+ xfin2)
        
        println("Calculated compositions for $(length(EA1)) tie lines")
        
        # Convert to cartesian coordinates (handle potential issues)
        x1 = Float64[]
        y1 = Float64[]
        x2 = Float64[]
        y2 = Float64[]
        try
            for i in 1:length(WA1)
                if !any(isnan.([WA1[i], EA1[i], FUR1[i]])) && !any(isinf.([WA1[i], EA1[i], FUR1[i]]))
                    xi, yi = ternary_to_cartesian(WA1[i], EA1[i], FUR1[i])
                    push!(x1, xi)
                    push!(y1, yi)
                end
            end
            for i in 1:length(WA2)
                if !any(isnan.([WA2[i], EA2[i], FUR2[i]])) && !any(isinf.([WA2[i], EA2[i], FUR2[i]]))
                    xi, yi = ternary_to_cartesian(WA2[i], EA2[i], FUR2[i])
                    push!(x2, xi)
                    push!(y2, yi)
                end
            end
            println("Converted $(length(x1)) EA-rich points and $(length(x2)) water-rich points to cartesian")
            if length(x1) == 0 || length(x2) == 0
                println("Warning: No valid points after conversion - skipping this file")
                continue
            end
        catch e
            println("Error in coordinate conversion: $e")
            continue
        end

        # Select 6 representative tie lines (evenly spaced) - use valid points only
        n_points = min(length(x1), length(x2))
        n_reps = min(6, n_points)
        if n_points == 0
            indices = Int[]
        else
            # Use floor to ensure indices are unique and within bounds, and always include first and last
            raw_indices = floor.(Int, range(1, n_points, length=n_reps))
            indices = unique(clamp.(raw_indices, 1, n_points))
        end
        
        # Create the plot
        plt = plot(size=(800, 700), dpi=300)
        
        # Draw triangle and grid (no ticks, clean axes)
        draw_ternary_triangle_clean!()
        add_ternary_grid!(0.1)
        
        # Plot all phase points as small dots
        scatter!(x1, y1,
                markersize=2,
                markercolor=:lightblue,
                markerstrokewidth=0,
                alpha=0.7,
                label="EA-rich phase (all points)")
        
        scatter!(x2, y2,
                markersize=2,
                markercolor=:lightcoral,
                markerstrokewidth=0,
                alpha=0.7,
                label="Water-rich phase (all points)")
        
        # Plot representative tie lines
        for (i, idx) in enumerate(indices)
            plot!([x1[idx], x2[idx]], [y1[idx], y2[idx]],
                  color=:blue,
                  linestyle=:dash,
                  linewidth=2,
                  label=(i==1 ? "Representative tie lines" : ""))
        end
        
        # Plot selected tie line endpoints with larger markers
        scatter!(x1[indices], y1[indices],
                markershape=:square,
                markersize=8,
                markercolor=:white,
                markerstrokecolor=:black,
                markerstrokewidth=2,
                label="EA-rich (selected)")
        
        scatter!(x2[indices], y2[indices],
                markershape=:square,
                markersize=8,
                markercolor=:black,
                markerstrokecolor=:black,
                markerstrokewidth=2,
                label="Water-rich (selected)")
        
        # Add title and format
        plot!(title="Water-Ethyl Acetate-Furfural Ternary Phase Diagram\\n$csv_file ($n_points tie lines, $(length(indices)) shown)",
              titlefontsize=12,
              legend=:topright,
              legendfontsize=8)
        

        
        # Print some statistics
        println("\\n=== PHASE DIAGRAM STATISTICS ===")
        println("Total tie lines: $(length(EA1))")
        println("Selected tie lines at indices: $indices")
        println("\\nEA-rich phase composition ranges:")
        println("  Water: $(round(minimum(WA1), digits=3)) - $(round(maximum(WA1), digits=3))")
        println("  EA: $(round(minimum(EA1), digits=3)) - $(round(maximum(EA1), digits=3))")
        println("  Furfural: $(round(minimum(FUR1), digits=3)) - $(round(maximum(FUR1), digits=3))")
        println("\\nWater-rich phase composition ranges:")
        println("  Water: $(round(minimum(WA2), digits=3)) - $(round(maximum(WA2), digits=3))")
        println("  EA: $(round(minimum(EA2), digits=3)) - $(round(maximum(EA2), digits=3))")
        println("  Furfural: $(round(minimum(FUR2), digits=3)) - $(round(maximum(FUR2), digits=3))")
        
        # Verify compositions sum to 1 (quality check)
        println("\\n=== DATA QUALITY CHECK ===")
        ea1_sums = EA1 .+ WA1 .+ FUR1
        ea2_sums = EA2 .+ WA2 .+ FUR2
        println("EA-rich phase sum check - min: $(round(minimum(ea1_sums), digits=4)), max: $(round(maximum(ea1_sums), digits=4))")
        println("Water-rich phase sum check - min: $(round(minimum(ea2_sums), digits=4)), max: $(round(maximum(ea2_sums), digits=4))")
        
    catch e
        println("Error loading or processing data from $csv_file: ", e)
        println("\\nPlease check:")
        println("1. CSV file contains required columns: XEAIN1, XWIN1, XFIN1, XEAIN2, XWIN2, XFIN2")
        println("2. Data format (numeric values, no missing data)")
        continue
    end
end

println("\\n" * "="^50)
println("Processing completed!")
println("="^50)