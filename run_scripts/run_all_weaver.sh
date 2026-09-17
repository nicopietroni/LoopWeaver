#!/bin/bash

set -m   # enable job control

cleanup() {
    echo "🛑 Batch interrupted, killing all child processes"
    kill -- -$$ 2>/dev/null
    exit 130
}
trap cleanup INT TERM

# ─────────────────────────────────────────────
# CONFIGURATION
#   Optional argument: dataset base folder containing the four variant
#   subfolders. Defaults to the full Flowloops dataset.
# ─────────────────────────────────────────────
DATASET_BASE="${1:-/Users/edwardsu/Projects/phd/Flowloops_dataset}"

if [ ! -d "$DATASET_BASE" ]; then
    echo "Error: dataset base folder not found: '$DATASET_BASE'"
    echo "Usage: $0 [dataset_base_folder]"
    exit 1
fi
DATASET_BASE="$(cd "$DATASET_BASE" && pwd)"

DIR_FEATURES_NO_SYMM="$DATASET_BASE/no_symm_features"
DIR_FEATURES_SYMM="$DATASET_BASE/symm_features"
DIR_NO_FEATURES_NO_SYMM="$DATASET_BASE/no_symm_no_features"
DIR_NO_FEATURES_SYMM="$DATASET_BASE/symm_no_features"

# ERROR_SETUP_VALUES=(0 1 2 3)
ERROR_SETUP_VALUES=(1)

# Locate the loop_weaver executable relative to this script
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_DIR="$SCRIPT_DIR/../build"
WEAVER="$BUILD_DIR/loop_weaver"

if [ ! -x "$WEAVER" ]; then
    echo "Error: loop_weaver executable not found at $WEAVER (build it first)"
    exit 1
fi

# Work from build so all outputs are created there
cd "$BUILD_DIR"

TIMEOUT=1200
CALL_DIR="$(pwd)"
OUTPUT_DIR="$CALL_DIR/output"

# ─────────────────────────────────────────────
# CORE FUNCTION
#   run_variant <folder> <angle> <use_symm> <variant_name> <error_setup>
#
#   Output structure:
#     <OUTPUT_DIR>/<variant_name>/es<error_setup>/  <- output files land here
#     <OUTPUT_DIR>/<variant_name>/timeouts.txt      <- one log per variant
# ─────────────────────────────────────────────
run_variant() {
    local folder="$1"
    local angle="$2"
    local use_symm="$3"    # "1" = add -symm, "0" = omit
    local variant_name="$4"
    local error_setup="$5"

    local TIMEOUT_LOG="${OUTPUT_DIR}/${variant_name}/timeouts.txt"
    local DONE_DIR="${OUTPUT_DIR}/${variant_name}/es${error_setup}"
    mkdir -p "$DONE_DIR"

    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "▶ Variant    : $variant_name"
    echo "  error_setup: $error_setup"
    echo "  Folder     : $folder"
    echo "  Output dir : $DONE_DIR"
    echo "  Started    : $(date)"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "=== es=${error_setup} run started at $(date) ===" >> "$TIMEOUT_LOG"

    shopt -s nullglob
    for file in "$folder"/*.obj; do
        echo "  Processing: $file"

        # Build the command arguments
        local extra_args=""
        [ "$use_symm" = "1" ] && extra_args="-symm"

        local start_time end_time elapsed status
        start_time=$(date +%s)

        gtimeout --kill-after=10s "$TIMEOUT" \
            "$WEAVER" "$file" \
            -error_setup "$error_setup" \
            -batch \
            -angle "$angle" \
            $extra_args

        status=$?
        end_time=$(date +%s)
        elapsed=$((end_time - start_time))

        if [ $status -eq 124 ]; then
            echo "  ⏱️  Timed out: $file (${elapsed}s)"
            echo "$(date '+%Y-%m-%d %H:%M:%S')  TIMEOUT  $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
            continue
        elif [ $status -ne 0 ]; then
            echo "  ❌ Error ($status): $file (${elapsed}s)"
            echo "$(date '+%Y-%m-%d %H:%M:%S')  ERROR    $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
            continue
        fi

        # Move output files whose name contains the input stem
        local stem moved out
        stem=$(basename "$file" .obj)
        moved=0
        for out in "$CALL_DIR"/*"${stem}"*; do
            [ -d "$out" ] && continue
            [[ "$out" == *timeouts* ]] && continue
            mv "$out" "$DONE_DIR/"
            moved=$((moved + 1))
        done

        echo "  ✅ Success: $file (${elapsed}s) → moved ${moved} output files to $DONE_DIR"
        echo "$(date '+%Y-%m-%d %H:%M:%S')  SUCCESS  $file  (${elapsed}s)" >> "$TIMEOUT_LOG"
    done

    echo "  Finished: $variant_name / es=${error_setup} at $(date)"
}

# ─────────────────────────────────────────────
# VALIDATE FOLDERS
# ─────────────────────────────────────────────
for dir_var in DIR_FEATURES_NO_SYMM DIR_FEATURES_SYMM DIR_NO_FEATURES_NO_SYMM DIR_NO_FEATURES_SYMM; do
    dir_path="${!dir_var}"
    if [ ! -d "$dir_path" ]; then
        echo "❌ Directory not found for $dir_var: '$dir_path'"
        echo "   Edit the CONFIGURATION section at the top of this script."
        exit 1
    fi
done

# ─────────────────────────────────────────────
# MAIN LOOP
#   Iterate through all error_setup values, running all 4 variants each time.
# ─────────────────────────────────────────────
echo "🚀 Starting all weaver variants — $(date)"
echo "   error_setup sequence: ${ERROR_SETUP_VALUES[*]}"

for es in "${ERROR_SETUP_VALUES[@]}"; do
    run_variant "$DIR_FEATURES_NO_SYMM"    45  0 "features_no_symm"    "$es"
    run_variant "$DIR_FEATURES_SYMM"       45  1 "features_symm"       "$es"
    run_variant "$DIR_NO_FEATURES_NO_SYMM" 180 0 "no_features_no_symm" "$es"
    run_variant "$DIR_NO_FEATURES_SYMM"    180 1 "no_features_symm"    "$es"
done

echo ""
echo "🏁 All variants and ccability values complete — $(date)"