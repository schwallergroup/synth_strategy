"""
Visualise Retrieval Results Script (Minimal Version)

This script generates images for synthesis routes based on the output of a retrieval
evaluation. It uses the `aizynthfinder` library for a robust and minimal
implementation of route image generation.

**Prerequisites:**
You must install `aizynthfinder` for this script to work. It will handle all
the necessary dependencies for image creation.

    pip install aizynthfinder

**Usage:**
    python visualise_retrieval_results.py \
        --retrieval_file /path/to/your/retrieval_results_model_X.json \
        --benchmark_file /path/to/your/merged_bench.json \
        --output_dir /path/to/your/visualizations \
        --top_n 5
"""
import argparse
import json
from pathlib import Path
from typing import Dict

from tqdm import tqdm

# This script requires the `aizynthfinder` package.
try:
    from aizynthfinder.reactiontree import ReactionTree
except ImportError:
    print("Error: The 'aizynthfinder' package is not installed.")
    print("Please install it to run this script: pip install aizynthfinder")
    exit(1)


def generate_route_image(route_dict: Dict, output_path: Path):
    """
    Generates and saves a route image using the aizynthfinder library.
    Returns True on success, False on failure.
    """
    if not route_dict:
        print(f"  - Warning: Skipping image generation for empty route_dict -> {output_path}")
        return False

    try:
        # 1. Create the ReactionTree object from the dictionary data
        reaction_tree = ReactionTree.from_dict(route_dict)

        # 2. Generate the PIL image object
        image = reaction_tree.to_image()

        # 3. Ensure the output directory exists and save the image
        output_path.parent.mkdir(parents=True, exist_ok=True)
        image.save(output_path)

    except Exception as e:
        print(f"  - Error: Failed to create or save image for {output_path.name}. Skipping. Details: {e}")
        return False
    return True

def create_route_lookup(benchmark_path: Path) -> Dict[str, Dict]:
    """Loads benchmark data and creates a fast route_id -> route_data mapping."""
    print(f"Loading benchmark data from '{benchmark_path}' to create a route lookup table...")
    with benchmark_path.open('r', encoding='utf-8') as f:
        benchmark_data = json.load(f)

    route_lookup = {item["selected_route_id"]: item["route_data"] for item in benchmark_data}
    print(f"Created lookup for {len(route_lookup)} routes.")
    return route_lookup

def main():
    """Main execution function to generate retrieval visualizations."""
    parser = argparse.ArgumentParser(
        description="Generate images for retrieved synthesis routes using aizynthfinder.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--retrieval_file",
        type=Path,
        required=True,
        help="Path to the JSON file containing retrieval results from the evaluation script."
    )
    parser.add_argument(
        "--benchmark_file",
        type=Path,
        required=True,
        help="Path to the benchmark JSON file containing the full `route_data` for all routes."
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        required=True,
        help="Directory to save the generated images. A subdirectory will be created for each model and query."
    )
    parser.add_argument(
        "--top_n",
        type=int,
        default=5,
        help="Number of top retrieved results to visualize for each query."
    )
    parser.add_argument(
        "--limit_queries",
        type=int,
        default=None,
        help="Limit the visualization to the first N queries in the retrieval file (for quick tests)."
    )
    parser.add_argument(
        "--image_format",
        type=str,
        default="png",
        help="The output format for the generated route images (e.g., png, jpg). Determined by extension."
    )
    args = parser.parse_args()

    # --- Validation and Setup ---
    if not args.retrieval_file.is_file():
        print(f"Error: Retrieval file not found at '{args.retrieval_file}'")
        return
    if not args.benchmark_file.is_file():
        print(f"Error: Benchmark file not found at '{args.benchmark_file}'")
        return

    model_name = args.retrieval_file.stem.replace("retrieval_results_", "")
    run_output_dir = args.output_dir / model_name
    run_output_dir.mkdir(parents=True, exist_ok=True)

    print("\n--- Starting Retrieval Visualization (Minimal Version) ---")
    print(f"  Model/Run Name: {model_name}")
    print(f"  Output Dir:     {run_output_dir}")
    print(f"  Visualize Top:  {args.top_n} results")
    if args.limit_queries:
        print(f"  Query Limit:    {args.limit_queries}")
    print("----------------------------------------------------------\n")

    # --- Data Loading ---
    route_lookup = create_route_lookup(args.benchmark_file)
    with args.retrieval_file.open('r', encoding='utf-8') as f:
        retrieval_data = json.load(f)

    queries_to_process = list(retrieval_data.items())
    if args.limit_queries:
        queries_to_process = queries_to_process[:args.limit_queries]

    # --- Image Generation Loop ---
    print(f"Generating images for {len(queries_to_process)} queries...")
    for query_id, retrieved_ids in tqdm(queries_to_process, desc="Processing Queries"):
        query_output_dir = run_output_dir / query_id

        # 1. Generate image for the query route
        query_route_data = route_lookup.get(query_id)
        if query_route_data:
            query_image_path = query_output_dir / f"query_{query_id}.{args.image_format}"
            generate_route_image(query_route_data, query_image_path)
        else:
            print(f"Warning: Could not find route data for query ID '{query_id}'")

        # 2. Generate images for the top N retrieved routes
        for i, retrieved_id in enumerate(retrieved_ids[:args.top_n]):
            rank = i + 1
            retrieved_route_data = route_lookup.get(retrieved_id)
            if retrieved_route_data:
                is_hit = (retrieved_id == query_id)
                hit_str = "HIT" if is_hit else "MISS"
                
                image_filename = f"rank_{rank:02d}_{hit_str}_{retrieved_id}.{args.image_format}"
                retrieved_image_path = query_output_dir / image_filename
                generate_route_image(retrieved_route_data, retrieved_image_path)
            else:
                 print(f"Warning: Could not find route data for retrieved ID '{retrieved_id}' at rank {rank}")

    print("\n--- Visualization Complete ---")
    print(f"All images have been saved to: {run_output_dir.resolve()}")
    print("------------------------------")

if __name__ == "__main__":
    main()