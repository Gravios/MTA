#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 input_file html_tag"
    exit 1
fi

# Input file and tag
input_file="$1"
html_tag="$2"

# Output file with .txt extension
output_file="${input_file%.*}.txt"

# Use awk to extract content between specified tags
awk -v tag="$html_tag" '
    BEGIN { RS = "<"tag">|</"tag">"; FS = "" }
    NR % 2 == 0 { print $0 }
' "$input_file" > "$output_file"

# Add newline after each entry if not already present
sed -i -e '$a\' "$output_file"

echo "Contents extracted to $output_file"
