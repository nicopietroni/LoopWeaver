#!/bin/bash


# Check argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <folder>"
    exit 1
fi

FOLDER="$1"

# Check folder exists
if [ ! -d "$FOLDER" ]; then
    echo "Error: '$FOLDER' is not a directory"
    exit 1
fi

# Loop over all .obj files in the folder
shopt -s nullglob
for file in "$FOLDER"/*.obj; do
    echo "Processing: $file"
    ./loop_weaver "$file" -batch -symm -angle 45
done

echo "All done."

