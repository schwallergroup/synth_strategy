#!/bin/bash
#SBATCH --job-name=route_annotation
#SBATCH --output=slurm_logs/annotation_%A_%a.out
#SBATCH --error=slurm_logs/annotation_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=1-5000 #<-- IMPORTANT: Set this to the total number of files to process

#======================================================================================
#
# SLURM Job Submission Script for Parallel Route Annotation
#
# Purpose:
#   This script is designed to be run as a SLURM job array. It processes a
#   collection of retrosynthesis route files (e.g., from AiZynthFinder) in parallel.
#   Each task in the array processes one input file using 'annotation.py',
#   which applies custom strategic bond disconnection logic.
#
# Prerequisites:
#   - A Python environment (e.g., Conda) with all necessary packages installed.
#   - The 'annotation.py' script.
#   - A directory of custom strategy function Python files.
#   - A directory of raw route files, numbered sequentially.
#
# Usage:
#   1. Configure the variables in the "USER CONFIGURATION" section below.
#   2. Adjust the SBATCH resource directives (e.g., --mem, --time, --array) as needed.
#   3. Submit the script to SLURM: sbatch submit_annotation.sh
#
#======================================================================================

# --- USER CONFIGURATION ---
# Please modify the following paths and settings to match your environment.

# 1. Project Directory (base path for your project)
PROJECT_DIR="/path/to/your/project"

# 2. Input Directory: Contains the raw route files to be annotated.
#    Files are expected to be named sequentially (e.g., smiles0001.json, smiles0002.json)
RAW_ROUTES_DIR="${PROJECT_DIR}/data/raw_routes"

# 3. Output Directory: Where the annotated route files will be saved.
ANNOTATED_ROUTES_DIR="${PROJECT_DIR}/data/annotated_routes"

# 4. Functions Directory: Contains the Python files with your custom strategy functions.
STRATEGY_FUNCTIONS_DIR="${PROJECT_DIR}/src/strategy_functions"

# 5. Conda Environment: Name of the Conda environment to activate.
CONDA_ENV_NAME="your-conda-env-name"

# --- END OF USER CONFIGURATION ---


# --- Script Setup ---
# This section sets up the environment and creates necessary directories.
# It is recommended not to modify this section.
set -e # Exit immediately if a command exits with a non-zero status.
set -u # Treat unset variables as an error when substituting.

echo "================================================="
echo "SLURM Job ID:        ${SLURM_JOB_ID}"
echo "SLURM Array Job ID:  ${SLURM_ARRAY_JOB_ID}"
echo "SLURM Array Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "================================================="

echo "Creating output directories if they do not exist..."
mkdir -p "${ANNOTATED_ROUTES_DIR}"
mkdir -p "slurm_logs" # Directory for SLURM's own log files


# --- Environment Activation ---
echo "Loading Conda environment: ${CONDA_ENV_NAME}"
# IMPORTANT: The path to 'conda.sh' may vary.
# Common locations: ~/miniconda3/etc/profile.d/conda.sh or ~/anaconda3/etc/profile.d/conda.sh
source "/path/to/your/miniconda3/etc/profile.d/conda.sh"
conda activate "${CONDA_ENV_NAME}"


# --- Task-specific Logic ---
# Format the task ID with leading zeros (e.g., 1 -> 0001). This helps match filenames.
# The '4' in '%04d' specifies a total width of 4 digits. Adjust if needed.
TASK_ID=$(printf "%04d" "${SLURM_ARRAY_TASK_ID}")
echo "Formatted Task ID: ${TASK_ID}"

# Define input and output file paths for this specific task
INPUT_FILE="${RAW_ROUTES_DIR}/smiles${TASK_ID}.json"
OUTPUT_FILE="${ANNOTATED_ROUTES_DIR}/annotated_smiles${TASK_ID}.json"

# Robustness Check: Ensure the input file for this task actually exists.
if [[ ! -f "${INPUT_FILE}" ]]; then
    echo "Error: Input file ${INPUT_FILE} not found! Skipping this task."
    exit 1
fi


# --- Execute Annotation Script ---
echo "--- Starting Route Annotation ---"
echo "Input File:      ${INPUT_FILE}"
echo "Output File:     ${OUTPUT_FILE}"
echo "Functions Dir:   ${STRATEGY_FUNCTIONS_DIR}"
echo "CPUs for task:   ${SLURM_CPUS_PER_TASK}"

python annotation.py \
    --input "${INPUT_FILE}" \
    --output "${OUTPUT_FILE}" \
    --functions-dir "${STRATEGY_FUNCTIONS_DIR}" \
    --workers "${SLURM_CPUS_PER_TASK}"

echo "--- Annotation for task ${TASK_ID} completed successfully. ---"