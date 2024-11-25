#!/bin/bash
calculate_percentiles() {
    local column=$1
    local file=$2
    local q1=$(awk -v col="$column" 'NR > 1 {print $col}' "$file" | sort -n | awk 'NR == int(NR * 0.25)')
    local q3=$(awk -v col="$column" 'NR > 1 {print $col}' "$file" | sort -n | awk 'NR == int(NR * 0.75)')
    echo "$q1 $q3"
}

# Input matrix file
input_file="Genepy_matrix.txt"
output_file="scaled_output.txt"

# Read the header and write it to the output file
head -n 1 "$input_file" > "$output_file"

# Process each column (starting from column 2)
num_columns=$(awk '{print NF}' "$input_file" | sort -nu | tail -n 1)
awk -v num_columns="$num_columns" 'NR > 1 {
    # First print the sample name (first column)
    printf "%s", $1;
    
    for (i=2; i<=num_columns; i++) {
        # Calculate Q1 and Q3 for each column
        cmd = "awk -v col="i" \x27NR > 1 {print $"col"}\x27 " input_matrix.txt " | sort -n | awk \x27NR == int(NR * 0.25)\x27"
        cmd | getline q1;
        close(cmd);
        cmd = "awk -v col="i" \x27NR > 1 {print $"col"}\x27 " input_matrix.txt " | sort -n | awk \x27NR == int(NR * 0.75)\x27"
        cmd | getline q3;
        close(cmd);
        
        # Min-max scaling for the value in column 'i'
        scaled = ($i - q1) / (q3 - q1);
        if (scaled < 0) scaled = 0;  # Ensure scaled values don't go below 0
        if (scaled > 1) scaled = 1;  # Ensure scaled values don't exceed 1

        # Print the scaled value
        printf " %.4f", scaled;
    }

    print "";  # Newline after each row
}' "$input_file" >> "$output_file"
echo "Scaling complete. Scaled output saved to $output_file."
