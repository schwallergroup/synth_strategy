#!/bin/bash
#======================================================================================
#
# Local Execution Script for Route Annotation and Clustering
#
# Purpose:
#   This script orchestrates a two-step analysis pipeline WITHOUT using SLURM.
#   It is a direct replacement for 'submit_annotation.sh' and 'submit_clustering.sh'.
#
#   Step 1: Annotation
#     - It iterates through all raw route files.
#     - It runs the 'annotation.py' script for each file in parallel, up to a
#       user-defined limit (MAX_PARALLEL_JOBS).
#
#   Step 2: Clustering
#     - After all annotation jobs are complete, it runs 'clustering.py' once to
#       aggregate the results and produce a final analysis file.
#
# Usage:
#   1. Place this script, 'annotation.py', and 'clustering.py' in your project's
#      main directory or a 'scripts' subdirectory.
#   2. Configure ALL variables in the "USER CONFIGURATION" section below.
#   3. Make the script executable: chmod +x run_local.sh
#   4. Run the script from your terminal: ./run_local.sh
#
#======================================================================================

# --- USER CONFIGURATION ---
# Please modify the following paths and settings to match your environment.

# 1. Project Directory (base path for your project)
PROJECT_DIR="/path/to/your/project"

# 2. Total Number of Files to Process
#    This replaces the '--array=1-5000' SLURM directive.
TOTAL_FILES=5000

# 3. Parallelism Settings for the Annotation Step
#    How many annotation processes to run at the same time.
#    Set this based on your machine's core count. A good starting point is half the number of your CPU cores.
MAX_PARALLEL_JOBS=8
#    How many CPUs each individual annotation process can use for its internal multiprocessing.
CPUS_PER_TASK=4 # This was '--cpus-per-task' in the SLURM script.

# 4. Input/Output Directories (same as in the original scripts)
RAW_ROUTES_DIR="${PROJECT_DIR}/data/raw_routes"
ANNOTATED_ROUTES_DIR="${PROJECT_DIR}/data/annotated_routes"
RESULTS_DIR="${PROJECT_DIR}/results"
STRATEGY_FUNCTIONS_DIR="${PROJECT_DIR}/src/strategy_functions"

# 5. Output Filename for the final analysis
FINAL_ANALYSIS_FILE="${RESULTS_DIR}/strategy_clustering_analysis.json"

# 6. Conda Environment
CONDA_ENV_NAME="your-conda-env-name"

# --- END OF USER CONFIGURATION ---


# --- Script Setup ---
set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error when substituting.

echo "================================================="
echo "Starting Local Annotation & Clustering Pipeline"
echo "Project Directory: ${PROJECT_DIR}"
echo "Total files to process: ${TOTAL_FILES}"
echo "Max parallel annotation jobs: ${MAX_PARALLEL_JOBS}"
echo "================================================="

echo "Creating output directories if they do not exist..."
mkdir -p "${ANNOTATED_ROUTES_DIR}"
mkdir -p "${RESULTS_DIR}"
mkdir -p "logs" # Create a local logs directory


# --- Environment Activation ---
echo "Loading Conda environment: ${CONDA_ENV_NAME}"
# IMPORTANT: The path to 'conda.sh' may vary.
# Common locations: ~/miniconda3/etc/profile.d/conda.sh or ~/anaconda3/etc/profile.d/conda.sh
source "/path/to/your/miniconda3/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"
echo "Conda environment activated."


# --- PART 1: PARALLEL ANNOTATION ---
echo
echo "--- Starting Part 1: Route Annotation (in parallel) ---"

# This loop replaces the SLURM job array.
# It runs up to MAX_PARALLEL_JOBS processes in the background.
for i in $(seq 1 ${TOTAL_FILES}); do
    # This block runs in the background for each file
    (
        # Format the task ID with leading zeros (e.g., 1 -> 0001)
        TASK_ID=$(printf "%04d" "${i}")
        
        # Define input and output file paths for this specific task
        INPUT_FILE="${RAW_ROUTES_DIR}/smiles${TASK_ID}.json"
        OUTPUT_FILE="${ANNOTATED_ROUTES_DIR}/annotated_smiles${TASK_ID}.json"
        LOG_FILE="logs/annotation_${TASK_ID}.log"

        echo "[Task ${TASK_ID}] Starting..."

        # Check if the input file exists before processing
        if [[ ! -f "${INPUT_FILE}" ]]; then
            echo "[Task ${TASK_ID}] Error: Input file ${INPUT_FILE} not found! Skipping."
            exit 1 # Exit the subshell, not the main script
        fi

        # Execute the annotation script, redirecting its output to a log file
        python annotation.py \
            --input "${INPUT_FILE}" \
            --output "${OUTPUT_FILE}" \
            --functions-dir "${STRATEGY_FUNCTIONS_DIR}" \
            --workers "${CPUS_PER_TASK}" > "${LOG_FILE}" 2>&1

        echo "[Task ${TASK_ID}] Completed successfully. Log: ${LOG_FILE}"

    ) & # The '&' runs the subshell in the background

    # Limit the number of concurrent jobs
    # When the number of background jobs reaches our limit, wait for one to finish
    if [[ $(jobs -r -p | wc -l) -ge ${MAX_PARALLEL_JOBS} ]]; then
        wait -n # Waits for the next job to terminate
    fi
done

# Wait for all remaining background jobs to finish before proceeding
echo "Waiting for all annotation tasks to complete..."
wait
echo "--- Part 1: Annotation complete. ---"


# --- PART 2: STRATEGY CLUSTERING ---
echo
echo "--- Starting Part 2: Strategy Clustering Analysis ---"
echo "Input Directory:  ${ANNOTATED_ROUTES_DIR}"
echo "Output File:      ${FINAL_ANALYSIS_FILE}"
echo "Code Source Dir:  ${STRATEGY_FUNCTIONS_DIR}"

python clustering.py \
    --input-dir "${ANNOTATED_ROUTES_DIR}" \
    --output "${FINAL_ANALYSIS_FILE}" \
    --code-dir "${STRATEGY_FUNCTIONS_DIR}"

echo "--- Part 2: Clustering finished successfully. ---"
echo "================================================="
echo "Pipeline complete!"
echo "Final results saved to: ${FINAL_ANALYSIS_FILE}"
echo "================================================="