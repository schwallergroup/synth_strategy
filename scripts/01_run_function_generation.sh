#!/bin/bash

# =================================================================
# Step 1: Strategy Function Generation (SynthStrategy)
# =================================================================
# Runs the LLM-based strategy extraction and code generation step
# over a batch of synthesis routes. Edit the configuration below to
# point at your own route file and to scale up N_SAMPLES for a full
# run; the defaults here are a small batch suitable for a quick
# end-to-end check of the pipeline.
#
# Configuration:
# - Input: first N_SAMPLES routes from INPUT_FILE.
# - Output: a JSON file in OUTPUT_DIR (see scripts/readme.md Step 1).
# - Model: MODEL_NAME (defaults to a fast model; use a stronger model
#   for production runs).
# =================================================================

# --- Configuration ---
# Path to the main Python script to be executed.
MAIN_SCRIPT="../src/synth_strategy/llm/function_generation.py"

# Input and Output paths (relative to the `scripts` directory).
INPUT_FILE="../data/routes/ref_routes_n1.json"
OUTPUT_DIR="../data/function_lib_test"

# Paths to the required chemical pattern dictionaries.
# Adjust these if your dictionary files are in a different location.
DICTIONARY_DIR="../data/patterns"
FG_PATH="${DICTIONARY_DIR}/functional_groups.json"
REACTION_PATH="${DICTIONARY_DIR}/smirks.json"
RING_PATH="${DICTIONARY_DIR}/chemical_rings_smiles.json"

# Processing parameters.
N_SAMPLES=5
START_INDEX=0
MODEL_NAME="google/gemini-2.5-flash" # A fast and capable model for testing.
MAX_CONCURRENT=5 # Lower concurrency for a small test to avoid rate limits.

# --- Script Execution ---

# Exit immediately if a command exits with a non-zero status.
set -e

echo "--- Starting Function Generation Run ---"

# Create the output directory if it doesn't exist.
# The -p flag ensures no error is thrown if the directory already exists.
echo "Ensuring output directory exists: ${OUTPUT_DIR}"
mkdir -p "$OUTPUT_DIR"

# Check if the main script exists
if [ ! -f "$MAIN_SCRIPT" ]; then
    echo "Error: Main script '$MAIN_SCRIPT' not found. Make sure you are in the correct directory."
    exit 1
fi

echo "Running the code generation script with the following parameters:"
echo "  - Input File:     ${INPUT_FILE}"
echo "  - Output Dir:     ${OUTPUT_DIR}"
echo "  - Model:          ${MODEL_NAME}"
echo "  - Samples:        ${N_SAMPLES} (starting from index ${START_INDEX})"
echo "---------------------------"

# Execute the main Python script with all arguments.
# Using backslashes for readability with multiple arguments.
python "$MAIN_SCRIPT" \
    --input-file "$INPUT_FILE" \
    --output-dir "$OUTPUT_DIR" \
    --model-name "$MODEL_NAME" \
    --n-samples "$N_SAMPLES" \
    --start-index "$START_INDEX" \
    --max-concurrent "$MAX_CONCURRENT" \
    --fg-path "$FG_PATH" \
    --reaction-path "$REACTION_PATH" \
    --ring-path "$RING_PATH"

echo "--- Function Generation Run Complete ---"
echo "Results have been saved in: ${OUTPUT_DIR}"
echo "You can check the output file named something like: train_t_${MODEL_NAME//\//_}_results.json"