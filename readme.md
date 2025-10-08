# Synthetic Strategy Function Pipeline

This repository contains a multi-stage pipeline for automatically generating, filtering, and refining Python functions that encode specific retrosynthesis strategies. The ultimate goal is to create a high-quality, structured library of these functions for use in route annotation, analysis, and steerable retrosynthesis models.

The workflow is orchestrated through a series of scripts located in the `scripts/` directory.

## Table of Contents

- [Workflow Overview](#workflow-overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [CLI Usage](#cli-usage)
  - [Available Commands](#available-commands)
  - [Command: annotate](#command-annotate)
  - [Command: cluster](#command-cluster)
  - [Command: retrieve](#command-retrieve)
  - [Command: visualize](#command-visualize)
  - [Common Workflows](#common-workflows)
  - [Configuration File](#configuration-file)
  - [Input Format](#input-format)
  - [Troubleshooting](#troubleshooting)
- [Step-by-Step Guide](#step-by-step-guide)
  - [Step 1: Function Generation](#step-1-function-generation)
  - [Step 2: Extract Functions to Python Files](#step-2-extract-functions-to-python-files)
  - [Step 3: Run the Filtering & Refinement Pipeline](#step-3-run-the-filtering--refinement-pipeline)
- [Advanced Usage: SLURM Job Submission](#advanced-usage-slurm-job-submission)
  - [Route Annotation](#route-annotation)
  - [Strategy Clustering](#strategy-clustering)
- [Local Execution (Non-SLURM)](#local-execution-non-slurm)
- [Directory Structure](#directory-structure)

## Workflow Overview

The entire process is designed to be run sequentially, where the output of one step becomes the input for the next.

### Function Library Creation (Scripts)

1. **Generate Functions** (`01_run_function_generation.sh`): An LLM is prompted with chemical reaction data to generate raw Python functions that attempt to identify synthetic strategies.

2. **Extract Functions** (`02_extract_functions_from_json.py`): The raw JSON output from the generation step is parsed, and each valid function is saved into its own `.py` file.

3. **Filter & Refine** (`03_run_filtering_pipeline.py`): This is the core multi-stage pipeline that uses LLMs to iteratively improve the quality of the generated functions:
   - **Stage 1 (Flash Filter)**: A fast, inexpensive LLM performs an initial quality check, removing obvious failures.
   - **Stage 2 (Pro Filter)**: A more powerful LLM performs a deeper analysis, fixing minor bugs and improving descriptions.
   - **Stage 3 (Metadata Generation)**: The logic of each passing function is deconstructed into a structured JSON format.
   - **Stage 4 (Return Modification)**: The function's code is automatically refactored to return not just a boolean, but the structured metadata of what it found.

### Using the Function Library (CLI)

4. **Annotate Routes** (`synth-strategy annotate`): Use the CLI to apply the function library to synthesis routes, identifying which strategies are present in each route.

5. **Cluster Strategies** (`synth-strategy cluster`): Perform clustering analysis to discover common strategy patterns and their relationships.

6. **Visualize Results** (`synth-strategy visualize`): Generate visualizations to understand strategy distributions and clustering results.

> **Note**: The CLI (`synth-strategy`) is the recommended way to use the function library for annotation, clustering, and visualization. For large-scale HPC processing, SLURM scripts are also available (see [Advanced Usage](#advanced-usage-slurm-job-submission)).

## Prerequisites

- Python 3.9+
- A Conda or other virtual environment management tool.
- An API key for an LLM provider (e.g., OpenRouter) set as an environment variable:


## Installation

### 1. Clone the Repository

```bash
git clone https://github.com/yourusername/synth-strategy-pipeline.git
cd synth-strategy-pipeline
```

### 2. Set Up Environment

Create and activate a conda environment:

```bash
conda create --name synth_strategy_env python=3.10
conda activate synth-strategy
conda install -c conda-forge rdkit
```

### 3. Install Dependencies

Install as required:

```bash
pip install --upgrade -e .
```

### 4. Configure API Keys

The pipeline requires API keys for LLM providers. Set these as environment variables:

```bash
# OpenRouter API key (required for function generation and filtering)
export OPENROUTER_API_KEY="your_openrouter_api_key_here"

```

To make these permanent, add them to your shell configuration file (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
echo 'export OPENROUTER_API_KEY="your_openrouter_api_key_here"' >> ~/.bashrc
source ~/.bashrc
```

### 5. Download Pre-annotated Data (Optional)

To get started quickly with pre-annotated USPTO route data, run the download script:

```bash
python scripts/download_uspto_data.py
```

This will:
- Download 4 annotated route files from Figshare (~2.5GB total)
- Extract them to `data/uspto_st/`
- Verify file integrity

The downloaded data includes thousands of annotated synthesis routes ready for retrieval and clustering analysis.

**Manual Download:**
If you prefer to download manually, the data is available at:
https://figshare.com/account/articles/30146374

## CLI Usage

The `synth-strategy` CLI provides a unified interface for route annotation, clustering, retrieval, and visualization. After installation, the CLI is available globally as `synth-strategy`.

### Available Commands

- `annotate` - Annotate routes with strategy functions
- `cluster` - Perform strategy clustering analysis
- `retrieve` - Search for strategies using natural language queries
- `visualize` - Generate visualizations from results

### Command: annotate

Annotate synthesis routes by applying strategy functions to identify which strategies are present in each route.

**Basic Usage:**

```bash
# Annotate routes from a directory
synth-strategy annotate \
  --input-dir data/routes/ \
  --functions-dir data/strategy_function_library/ \
  --output results/

# Annotate routes from a single file
synth-strategy annotate \
  --input-file data/test/test_routes.json \
  --functions-dir data/strategy_function_library/ \
  --output results/
```

**Advanced Usage:**

```bash
# Annotate and cluster in one step
synth-strategy annotate \
  --input-dir data/routes/ \
  --functions-dir data/strategy_function_library/ \
  --output results/ \
  --cluster \
  --visualize \
  --interactive

# Use multiple workers for faster processing
synth-strategy annotate \
  --input-file data/routes.json \
  --functions-dir data/strategy_function_library/ \
  --workers 8 \
  --output results/
```

**Options:**

- `--input-dir` - Directory containing route JSON files
- `--input-file` - Single route JSON file
- `--functions-dir` - Directory containing strategy function Python files
- `--output` / `-o` - Output directory for results
- `--workers` - Number of parallel workers (default: CPU count)
- `--cluster` - Run clustering after annotation
- `--visualize` - Generate visualizations (requires --cluster)
- `--interactive` - Create interactive visualizations
- `--config` / `-c` - Configuration file path

**Output:**

Creates `annotated_routes.json` with routes annotated with:
- `passing_functions` - Strategy functions that matched the route
- `errored_functions` - Functions that encountered errors

### Command: cluster

Perform clustering analysis on annotated routes to identify common strategy patterns.

**Basic Usage:**

```bash
# Cluster annotated routes
synth-strategy cluster \
  --input-file results/annotated_routes.json \
  --functions-dir data/strategy_function_library/ \
  --output results/

# Cluster and visualize
synth-strategy cluster \
  --input-dir data/annotated_routes/ \
  --functions-dir data/strategy_function_library/ \
  --output results/ \
  --visualize
```

**Options:**

- `--input-dir` - Directory containing annotated route JSON files
- `--input-file` - Single annotated route JSON file
- `--functions-dir` - Directory containing strategy function source code (for docstring extraction)
- `--output` / `-o` - Output directory for results
- `--visualize` - Generate visualizations after clustering

**Output:**

Creates `clustering_results.json` containing:
- `optimal_k` - Optimal number of clusters found
- `cluster_allocations` - Mapping of routes to cluster IDs
- `cluster_defining_features` - Key strategies that define each cluster
- `function_docstrings` - Documentation for each strategy function

### Command: retrieve

Search for synthesis strategies using natural language queries or structured queries.

**Basic Usage:**

```bash
# Simple text query
synth-strategy retrieve \
  --query "oxidation strategy for alcohol to aldehyde" \
  --metadata-db data/metadata.db \
  --route-db data/routes/ \
  --top-k 10 \
  --output results/

# Complex query from JSON file
synth-strategy retrieve \
  --query-file queries/complex_query.json \
  --metadata-db data/metadata.db \
  --route-db data/routes/ \
  --top-k 20 \
  --output results/ \
  --visualize
```

**Options:**

- `--query` - Natural language text query
- `--query-file` - JSON file containing structured query
- `--metadata-db` - Path to metadata database
- `--route-db` - Path to route database directory
- `--embedding-cache` - Path to embedding cache (optional)
- `--top-k` - Number of results to return
- `--output` / `-o` - Output directory
- `--visualize` - Generate visualizations

**Output:**

Creates `retrieval_results.json` with ranked strategy matches.

#### Understanding the Retrieval System

The retrieval system uses a sophisticated two-stage process to find synthesis routes matching your query:

**Stage 1: Query Rewriting (LLM-powered)**

Your natural language query is automatically converted to a structured JSON format by an LLM (default: `google/gemini-2.5-pro`). This structured query contains:
- A natural language description for semantic search
- Exact filters for precise matching using controlled vocabularies

Example query: *"late stage amide coupling, preservation of a piperidine core and early stage ring formation"*

Gets rewritten to:
```json
{
  "operator": "AND",
  "queries": [
    {
      "query": {
        "natural_language_description": "Detects late-stage amide coupling strategy",
        "filters": {
          "OR": {
            "named_reactions": [
              "Acyl chloride with primary amine to amide",
              "Carboxylic acid with primary amine to amide",
              "amide_formation"
            ]
          }
        }
      }
    },
    {
      "query": {
        "natural_language_description": "Preservation of piperidine core",
        "filters": {"ring_systems": ["piperidine"]}
      }
    },
    {
      "query": {
        "natural_language_description": "Early-stage ring formation",
        "filters": {"named_reactions": ["ring_formation"]}
      }
    }
  ]
}
```

**Stage 2: Two-Phase Retrieval**

For each sub-query, the system performs:

1. **Semantic Search (Broad Funnel)**
   - Uses embeddings to find functions with similar descriptions
   - Configurable via `--top-n-functions` (e.g., top 40 most similar)
   - Creates a large candidate pool

2. **Exact Filtering (Narrowing Down)**
   - Performs strict exact-match filtering on atomic checks
   - Uses controlled vocabularies: `named_reactions`, `ring_systems`, `functional_groups`
   - Only keeps functions that exactly match filter criteria

3. **Route Matching**
   - Finds routes containing the filtered functions
   - Checks instance-level matching (specific reaction details)
   - Scores routes by semantic similarity

4. **Ranking**
   - Primary: Number of sub-queries matched
   - Secondary: Average cosine similarity score
   - Returns top-k results

**Configuring Retrieval Behavior:**

The system is highly configurable through multiple mechanisms:

1. **Configuration File** (`config.yaml`):
```yaml
retrieval:
  embedding_model: "all-MiniLM-L6-v2"
  top_k: 10
  top_n_functions: 40  # null = use all functions
  metadata_db: "data/function_metadata_database.json"
  route_db_dir: "data/routes"
  embedding_cache: "data/retrieval_embeddings"
```

2. **CLI Override**:
```bash
synth-strategy retrieve \
  --query "your query" \
  --metadata-db path/to/metadata.json \
  --route-db path/to/routes/ \
  --embedding-cache path/to/cache.pkl \
  --top-k 20
```

3. **Python API**:
```python
from synth_strategy.retrieval.retriever import StrategyRetriever
from synth_strategy.models import SentenceTransformerEmbedder

embedder = SentenceTransformerEmbedder(model_name="all-MiniLM-L6-v2")
retriever = StrategyRetriever(
    metadata_db_path="data/function_metadata_database.json",
    route_db_dir="data/routes",
    embedding_cache_path="data/retrieval_embeddings/cache.pkl",
    embedder=embedder
)

results = retriever.retrieve_complex(query_dict, top_k=10, top_n_functions=40)
```

**Using a New Route Database:**

The system is designed to easily work with new route databases:

1. **Prepare Your Routes**: Routes must be annotated with strategy functions first:
```bash
synth-strategy annotate \
  --input-dir new_routes/ \
  --functions-dir data/strategy_functions/ \
  --output annotated_routes/
```

2. **Point to New Database**:
```bash
synth-strategy retrieve \
  --query "your query" \
  --route-db annotated_routes/ \
  --metadata-db data/function_metadata_database.json
```

3. **Automatic Processing**:
   - The system scans all `*.json` files in the directory
   - Builds an inverted index (function → routes mapping)
   - Caches the index for faster subsequent loads
   - Auto-invalidates cache when source files change

**Data Requirements:**

The retrieval system requires:
- **Metadata Database**: `function_metadata_database.json` containing function descriptions and atomic checks
- **Embedding Cache**: Pre-computed embeddings for all function descriptions
- **Route Database**: Directory of annotated route JSON files with `passing_functions` field
- **Inverted Index**: Auto-generated and cached in `.retriever_cache/`

**Performance Optimization:**

- **Caching**: Inverted index and embeddings are cached for fast repeated queries
- **Filtering**: Only functions that pass at least once in your route database are loaded
- **Parallel Processing**: Embedding computation can be parallelized

### Command: visualize

Generate visualizations from clustering or retrieval results.

**Usage:**

```bash
# Visualize clustering results
synth-strategy visualize \
  --type clustering \
  --input results/clustering_results.json \
  --annotated_dir data/annotated_routes/ \
  --output visualizations/

# Visualize retrieval results
synth-strategy visualize \
  --type retrieval \
  --input results/retrieval_results.json \
  --output visualizations/
```

**Options:**

- `--type` - Type of visualization: `clustering` or `retrieval` (required)
- `--input` / `-i` - Input results file (required)
- `--output` / `-o` - Output directory
- `--annotated_dir` - Directory containing annotated routes (required for clustering visualization)

### Common Workflows

**Workflow 1: Quick Analysis**

```bash
# Annotate, cluster, and visualize in one command
synth-strategy annotate \
  --input-file data/routes.json \
  --functions-dir data/strategy_function_library/ \
  --output results/ \
  --cluster \
  --visualize \
  --interactive
```

**Workflow 2: Step-by-Step Analysis**

```bash
# Step 1: Annotate routes
synth-strategy annotate \
  --input-dir data/routes/ \
  --functions-dir data/strategy_function_library/ \
  --output results/

# Step 2: Cluster annotated routes
synth-strategy cluster \
  --input-file results/annotated_routes.json \
  --functions-dir data/strategy_function_library/ \
  --output results/

# Step 3: Visualize results
synth-strategy visualize \
  --type clustering \
  --input results/clustering_results.json \
  --annotated_dir data/annotated_routes/ \
  --output visualizations/
```

**Workflow 3: Large-Scale Processing**

```bash
# Process large datasets with parallel workers
synth-strategy annotate \
  --input-dir data/large_dataset/ \
  --functions-dir data/strategy_function_library/ \
  --workers 16 \
  --output results/

# Cluster with visualization
synth-strategy cluster \
  --input-file results/annotated_routes.json \
  --functions-dir data/strategy_function_library/ \
  --output results/ \
  --visualize
```

### Configuration File

You can use a configuration file (YAML format) to set default parameters:

```yaml
# config.yaml
defaults:
  functions_dir: "data/strategy_function_library"
  output_dir: "results"
  workers: 8

retrieval:
  metadata_db: "data/metadata.db"
  route_db_dir: "data/routes"
  embedding_cache: "data/embeddings.pkl"
  top_k: 10
```

Use with any command:

```bash
synth-strategy annotate --config config.yaml --input-file data/routes.json
```

### Input Format

**Route JSON Format:**

Routes should be in JSON format as a list of route objects:

```json
[
  {
    "type": "mol",
    "smiles": "CC(=O)O",
    "metadata": {
      "target_smiles": "CC(=O)O"
    },
    "children": [
      {
        "type": "reaction",
        "metadata": {
          "mapped_reaction_smiles": "CCO>>CC(=O)O"
        },
        "children": [...]
      }
    ]
  }
]
```

### Troubleshooting

**Issue: "Not enough routes for clustering"**
- Ensure you have at least 2 annotated routes
- Check that routes have `passing_functions` populated

**Issue: "No passing functions found"**
- Verify that the strategy functions are compatible with your route format
- Check function error logs in the annotated routes

**Issue: "Module not found" errors**
- Ensure the package is installed: `pip install -e .`
- Activate the correct conda environment

## Step-by-Step Guide

All commands should be run from the `scripts/` directory.

### Step 1: Function Generation

This step uses an LLM to generate the initial set of raw strategy functions.

**Script:** `01_run_function_generation.sh`

This is a bash script that calls the underlying Python generation logic. Before running, you can configure the parameters inside the script, such as the input data, number of samples, and which LLM to use.

**How to Run:**

```bash
cd scripts/
./01_run_function_generation.sh
```

**Output:** A JSON file (e.g., `train_t_google_gemini-1.5-flash_results.json`) will be created in the output directory specified within the script (default: `../data/function_lib_test`).

### Step 2: Extract Functions to Python Files

This script parses the JSON output from Step 1 and creates individual `.py` files for each successfully generated function.

**Script:** `02_extract_functions_from_json.py`

**How to Run:**

You need to provide the path to the JSON file generated in the previous step and specify an output directory for the `.py` files.

```bash
# Replace with the actual path to your generated JSON
JSON_FILE="../data/function_lib_test/train_t_google_gemini-1.5-flash_results.json"

# This directory will contain all the raw Python functions
OUTPUT_DIR="../data/generated_functions_raw"

python 02_extract_functions_from_json.py "$JSON_FILE" --output-dir "$OUTPUT_DIR"
```

**Output:** The `../data/generated_functions_raw` directory will be populated with many `.py` files. This directory will be the main input for the next step.

### Step 3: Run the Filtering & Refinement Pipeline

This is the main orchestration script that takes the raw functions and puts them through the four-stage refinement process.

**Script:** `03_run_filtering_pipeline.py`

**How to Run:**

The script is highly configurable via command-line arguments. The most important one is `--source-code-dir`, which must point to the directory created in Step 2.

```bash
# Run the full pipeline on the raw functions
python 03_run_filtering_pipeline.py \
    --source-code-dir ../data/generated_functions_raw
```

**Key Parameters:**

- `--source-code-dir`: (Required) Path to the input directory of `.py` function files.
- `--base-output-dir`: Where the intermediate and final code artifacts will be stored. Default: `./pipeline_output`.
- `--reports-dir`: Where JSON reports for each stage will be saved. Default: `../data/filtering_data`.
- `--run-stages`: Specify which stages to run (e.g., `--run-stages 1 2` to only run the filtering). Default: `1 2 3 4`.
- `--flash-filter-level`: Choose what to keep from the first stage. `perfect` (default) or `good_and_perfect`.

**Output:**

- `./pipeline_output/`: This directory will contain the outputs of each stage in numbered subdirectories. The final, high-quality, and refactored functions will be located in `./pipeline_output/4_modified_functions_output/`.
- `../data/filtering_data/`: This will contain JSON reports detailing the results of each stage (e.g., `flash_filter_report.json`, `modification_report.json`).

## Advanced Usage: SLURM Job Submission

For large-scale processing, the final function library can be used to annotate routes on a high-performance computing (HPC) cluster using SLURM.

### Route Annotation

**Script:** `submit_annotation_jobs.sh`

This script submits a SLURM job array to process thousands of route files in parallel. Each task in the array annotates one route file with the entire library of strategy functions.

**Configuration:**

Before submitting, you must edit the "USER CONFIGURATION" section inside the script to set the correct paths for your project directory, input/output folders, and Conda environment.

**How to Run:**

```bash
sbatch submit_annotation_jobs.sh
```

### Strategy Clustering

**Script:** `submit_clustering_job.sh`

After the annotation jobs are complete, this script runs a final analysis job. It aggregates all the annotated route data and performs a clustering analysis to identify which strategies are most common and which tend to co-occur.

**Configuration:**

Similarly, you must edit the "USER CONFIGURATION" section to match your environment paths.

**How to Run:**

```bash
sbatch submit_clustering_job.sh
```

## Local Execution (Non-SLURM)

For users without access to a SLURM-based HPC cluster, a local execution script is provided that can run the annotation and clustering steps on a single machine.

**Script:** `run_local.sh`

This script orchestrates a two-step analysis pipeline without using SLURM:

1. **Step 1 - Annotation:** Iterates through all raw route files and runs the `annotation.py` script for each file in parallel, up to a user-defined limit (`MAX_PARALLEL_JOBS`).

2. **Step 2 - Clustering:** After all annotation jobs complete, runs `clustering.py` once to aggregate the results and produce a final analysis file.

### Configuration

Before running, you must edit the "USER CONFIGURATION" section inside the script to set:

- `PROJECT_DIR`: Base path for your project
- `TOTAL_FILES`: Number of files to process (replaces SLURM's `--array=1-5000`)
- `MAX_PARALLEL_JOBS`: How many annotation processes to run simultaneously (recommended: half of your CPU cores)
- `CPUS_PER_TASK`: How many CPUs each individual annotation process can use for its internal multiprocessing
- Input/output directories:
  - `RAW_ROUTES_DIR`: Location of input route files
  - `ANNOTATED_ROUTES_DIR`: Where annotated files will be saved
  - `RESULTS_DIR`: Where final analysis will be saved
  - `STRATEGY_FUNCTIONS_DIR`: Location of strategy functions library
- `CONDA_ENV_NAME`: Your conda environment name
- Conda installation path (update the `source` path in the script)

### How to Run

```bash
# 1. Make the script executable
chmod +x run_local.sh

# 2. Run the script
./run_local.sh
```

### Output

- **Annotation logs:** Individual log files for each task will be saved in `logs/annotation_XXXX.log`
- **Annotated routes:** Saved in the `ANNOTATED_ROUTES_DIR` directory
- **Final analysis:** `strategy_clustering_analysis.json` in the `RESULTS_DIR`

### Performance Notes

- The `MAX_PARALLEL_JOBS` setting controls memory usage and CPU load. Start conservatively and adjust based on your system's performance.
- Processing time will depend on your machine's specifications and the complexity of the route files.
- Monitor system resources during execution to optimize the parallelism settings for your hardware.

## Directory Structure

A typical project layout after running the pipeline:

```
.
├── data/
│   ├── function_lib_test/         # Raw JSON output from Step 1
│   ├── generated_functions_raw/   # Individual .py files from Step 2
│   └── filtering_data/            # Reports from the main pipeline (Step 3)
│
├── pipeline_output/
│   ├── 1_flash_filter_output/
│   │   └── functions/             # Functions that passed Stage 1
│   ├── 2_pro_filter_output/
│   │   └── functions/             # Functions that passed Stage 2
│   ├── 3_metadata_output/
│   │   └── function_metadata_database.json
│   └── 4_modified_functions_output/ # FINAL, REFINED FUNCTION LIBRARY
│
├── scripts/
│   ├── 01_run_function_generation.sh
│   ├── 02_extract_functions_from_json.py
│   ├── 03_run_filtering_pipeline.py
│   ├── submit_annotation_jobs.sh
│   └── submit_clustering_job.sh
│
└── src/
    └── synth_strategy/
        └── llm/                   # Core Python modules for the pipeline
```