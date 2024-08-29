#!/bin/bash

# Function to convert flashcards
convert_flashcards() {
    input_file=$1
    output_file=$2

    # Read the input file into an array
    mapfile -t lines < "$input_file"

    # Initialize an empty result array
    result=()

    # Process the lines in groups of three
    for ((i=0; i<${#lines[@]}; i+=3)); do
        question=$(echo "${lines[i+2]}" | tr -d '()')
        more_info=$(echo "${lines[i+1]}" | tr -d '[]')
        answer="${lines[i]}"
        if [[ -n "$question" ]]; then
            result+=("${question//-/} [$more_info] --- ${answer//-/}")
        fi
    done

    # Write the result to the output file
    printf "%s\n" "${result[@]}" > "$output_file"
}

# Check if the correct number of arguments is provided
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 input_file output_file"
    exit 1
fi

# Call the function with the provided arguments
convert_flashcards "$1" "$2"
