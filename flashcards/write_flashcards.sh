#!/bin/bash

# Default flashcards file
FLASHCARDS_FILE="flashcards.txt"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Parse options
while getopts "f:" opt; do
  case $opt in
    f) FLASHCARDS_FILE=$OPTARG ;;
    \?) echo "Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

while true; do
  echo -e "${YELLOW}Enter the question:${NC}"
  read question

  echo -e "${YELLOW}Enter the answer:${NC}"
  read answer

  flashcard="$question---$answer"
  echo "$flashcard" >> $FLASHCARDS_FILE
  echo -e "${GREEN}Flashcard added!${NC}"

  echo -e "${YELLOW}Press Enter to add another flashcard, or 'q' to quit.${NC}"
  read -n 1 choice
  if [ "$choice" == "q" ]; then
    break
  fi
done
