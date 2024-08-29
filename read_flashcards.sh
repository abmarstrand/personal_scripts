#!/bin/bash

# Default flashcards file
FLASHCARDS_FILE="flashcards.txt"

# Default values
NUM_FLASHCARDS=10
SUBSET=""

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Parse options
while getopts "n:s:f:" opt; do
  case $opt in
    n) NUM_FLASHCARDS=$OPTARG ;;
    s) SUBSET=$OPTARG ;;
    f) FLASHCARDS_FILE=$OPTARG ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done


# Function to read flashcards
read_flashcards() {
  local file=$1
  local num=$2
  local subset=$3
  local flashcards

  if [ -n "$subset" ]; then
    # Extract the specified range of lines
    IFS=':' read -r start end <<< "$subset"
    flashcards=$(sed -n "${start},${end}p" "$file")
  else
    # Read the entire file
    flashcards=$(cat "$file")
  fi

  # Select a random set of flashcards
  echo "$flashcards" | shuf -n "$num"
}

# Save flashcards to a variable
FLASHCARDS=$(read_flashcards "$FLASHCARDS_FILE" "$NUM_FLASHCARDS" "$SUBSET")

# Function to format flashcards
format_flashcards() {
  local flashcards=$1
  local formatted_flashcards=""

  while IFS= read -r line; do
    # Replace multiple spaces with a single space
    line=$(echo "$line" | tr -s ' ')
    # Split the line into question and answer
    IFS='---' read -r question answer <<< "$line"
    # Trim leading and trailing spaces
    question=$(echo "$question" | xargs)
    answer=$(echo "$answer" | xargs)
    answer=$(echo "$answer" | sed 's/^-*//;s/-*$//')
    # Append to formatted flashcards
    formatted_flashcards+=" $question\n $answer\n\n"
  done <<< "$flashcards"

  echo -e "$formatted_flashcards"
}

# Format the flashcards
FORMATTED_FLASHCARDS=$(format_flashcards "$FLASHCARDS")

# Function to display flashcards in TUI
display_flashcards() {
  local flashcards=$1
  local IFS=$'\n'
  local flashcard_array=($flashcards)
  local total_flashcards=$((${#flashcard_array[@]} / 2))

  for ((i=0; i<${#flashcard_array[@]}; i+=2)); do
    question="${flashcard_array[i]}"
    answer="${flashcard_array[i+1]}"
    current_flashcard=$((i / 2 + 1))

    # Display question
    clear
    echo -e "${CYAN}╔════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}║${NC} ${CYAN}Flashcard $current_flashcard/$total_flashcards${NC}"
    echo -e "${CYAN}║${NC} ${RED}QUESTION:${NC} $question"
    echo -e "${CYAN}╚════════════════════════════════════════════════╝${NC}"
    read -p "Press Enter to see the answer..."

    # Display question and answer
    clear
    echo -e "${CYAN}╔════════════════════════════════════════════════╗${NC}"
    echo -e "${CYAN}║${NC} ${CYAN}Flashcard $current_flashcard/$total_flashcards${NC}"
    echo -e "${CYAN}║${NC} ${RED}QUESTION:${NC} $question"
    echo -e "${CYAN}║${NC} ${GREEN}ANSWER:${NC} $answer"
    echo -e "${CYAN}╚════════════════════════════════════════════════╝${NC}"
    read -p "Press Enter to continue..."
  done
}
# Display the flashcards in TUI
display_flashcards "$FORMATTED_FLASHCARDS"
