#!/bin/bash

set -m   # enable job control

# Kill entire process group on Ctrl+C or termination
cleanup() {
    echo "🛑 Batch interrupted, killing all child processes"
    kill -- -$$ 2>/dev/null
    exit 130
}
trap cleanup INT TERM

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

# Settings
TIMEOUT=1800              # seconds (30 minutes)
TIMEOUT_LOG="timeouts.txt"

echo "=== Run started at $(date) ===" >> "$TIMEOUT_LOG"

# Loop over all .obj files in the folder
shopt -s nullglob
for file in "$FOLDER"/*.obj; do
    echo "Processing: $file"

    start_time=$(date +%s)
    timeout --foreground --kill-after=10s "$TIMEOUT" \
        ./loop_weaver "$file" -batch -symm -angle 45

    status=$?
    end_time=$(date +%s)
    elapsed=$((end_time - start_time))

    if [ $status -eq 124 ]; then
        echo "⏱️  Timed out: $file (${elapsed}s)"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  TIMEOUT  $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
        continue
    elif [ $status -ne 0 ]; then
        echo "❌ Error ($status): $file (${elapsed}s)"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  ERROR    $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
        continue
    fi
    echo "✅ Success: $file (${elapsed}s)"
    echo "$(date '+%Y-%m-%d %H:%M:%S')  SUCCESS  $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
done

echo "All done."
