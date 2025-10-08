"""
Visualise Clustering Results Script

This script takes the JSON output from `clustering.py` and generates a structured,
visual representation of the clusters. For each cluster, it creates a directory
containing:
1. A `summary.md` file describing the cluster's most defining strategic functions.
2. Image files for a sample of routes belonging to that cluster.

**Prerequisites:**
You must install `aizynthfinder` for this script to work.

    pip install aizynthfinder

**Usage:**
    python visualise_clustering_results.py \
        --analysis_file /path/to/strategy_clustering_analysis.json \
        --annotated_dir /path/to/annotated_routes/ \
        --output_dir /path/to/clustering_visuals/ \
        --examples_per_cluster 10
"""
import argparse
import json
import glob
from pathlib import Path
from typing import Dict, List, Any

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
        reaction_tree = ReactionTree.from_dict(route_dict)
        image = reaction_tree.to_image()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        image.save(output_path)
    except Exception as e:
        print(f"  - Error: Failed to create or save image for {output_path.name}. Skipping. Details: {e}")
        return False
    return True

def load_all_annotated_routes(annotated_dir: Path) -> List[Dict[str, Any]]:
    """
    Aggregates all routes from the JSON files in the specified directory.
    This is crucial as the indices in the analysis file refer to this aggregated list.
    """
    all_routes = []
    json_files = sorted(glob.glob(str(annotated_dir / '*.json'))) # Sorting is important for consistency
    if not json_files:
        raise FileNotFoundError(f"No JSON files found in the directory: {annotated_dir}")

    print(f"Loading and aggregating routes from {len(json_files)} files in '{annotated_dir}'...")
    for file_path in tqdm(json_files, desc="Loading annotated files"):
        try:
            with open(file_path, 'r') as f:
                all_routes.extend(json.load(f))
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load or parse {file_path}: {e}")
    
    print(f"Successfully loaded a total of {len(all_routes)} routes.")
    return all_routes

def create_cluster_summary(
    cluster_dir: Path,
    cluster_id: int,
    defining_features: List[List[Any]],
    docstrings: Dict[str, str],
    top_n_features: int = 5
):
    """Creates a Markdown summary file for a single cluster."""
    summary_path = cluster_dir / "summary.md"
    with summary_path.open('w') as f:
        f.write(f"# Cluster {cluster_id} Summary\n\n")
        f.write("This cluster is primarily defined by the following strategic functions, "
                "ordered by how distinctive they are to this group.\n\n")
        
        f.write(f"## Top {top_n_features} Defining Functions\n\n")

        # feature is a dict with"function", "score" and "description"
        # rewrite accordingly
        for feature in defining_features[:top_n_features]:
            feature_name = feature["function"]
            distinctiveness_score = feature["score"]
            docstring = feature["description"]
            f.write(f"### 1. `{feature_name}`\n")
            f.write(f"*Distinctiveness Score: {distinctiveness_score:.4f}*\n\n")
            f.write("```\n")
            f.write(docstring.strip() + "\n")
            f.write("```\n\n")
            f.write("---\n\n")

        

def main():
    parser = argparse.ArgumentParser(
        description="Generate visual summaries for route clustering results.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--analysis_file",
        type=Path,
        required=True,
        help="Path to the strategy_clustering_analysis.json file from clustering.py."
    )
    parser.add_argument(
        "--annotated_dir",
        type=Path,
        required=True,
        help="Directory containing the original ANNOTATED JSON files that were used as input for clustering."
    )
    parser.add_argument(
        "--output_dir",
        type=Path,
        required=True,
        help="Directory to save the generated cluster summaries and images."
    )
    parser.add_argument(
        "--examples_per_cluster",
        type=int,
        default=10,
        help="Number of example route images to generate for each cluster."
    )
    parser.add_argument(
        "--image_format",
        type=str,
        default="png",
        help="The output format for the generated route images (e.g., png, jpg)."
    )
    args = parser.parse_args()

    # --- Validation and Setup ---
    if not args.analysis_file.is_file():
        print(f"Error: Analysis file not found at '{args.analysis_file}'")
        return

    args.output_dir.mkdir(parents=True, exist_ok=True)
    
    print("\n--- Starting Clustering Visualization ---")
    print(f"  Analysis File:  {args.analysis_file}")
    print(f"  Annotated Dir:  {args.annotated_dir}")
    print(f"  Output Dir:     {args.output_dir}")
    print("-----------------------------------------\n")
    
    # --- Data Loading ---
    print(f"Loading analysis data from '{args.analysis_file}'...")
    with args.analysis_file.open('r') as f:
        analysis_data = json.load(f)

    if analysis_data.get('status') != 'success':
        print(f"Clustering status is '{analysis_data.get('status')}'. Reason: {analysis_data.get('reason')}. Exiting.")
        return

    all_routes = load_all_annotated_routes(args.annotated_dir)
    if not all_routes:
        print("Error: No routes were loaded from the annotated directory. Cannot proceed.")
        return
    
    # --- Generate Visualizations ---
    num_clusters = analysis_data['optimal_k']
    allocations = analysis_data['cluster_allocations']
    
    print(f"\nFound {num_clusters} clusters. Generating summaries and example images...")
    
    for cluster_id in range(num_clusters):
        cluster_dir = args.output_dir / f"cluster_{cluster_id}"
        cluster_dir.mkdir(exist_ok=True)

        print(f"\nProcessing Cluster {cluster_id}...")
        
        # 1. Create a human-readable summary file for the cluster
        create_cluster_summary(
            cluster_dir=cluster_dir,
            cluster_id=cluster_id,
            defining_features=analysis_data['cluster_defining_features'][str(cluster_id)],
            docstrings=analysis_data['function_docstrings']
        )
        
        # 2. Find all routes belonging to this cluster
        route_indices = [
            int(idx) for idx, c_id in allocations.items() if c_id == cluster_id
        ]
        
        print(f"  Found {len(route_indices)} routes in this cluster. Generating images for top {args.examples_per_cluster} examples.")

        # 3. Generate images for a sample of these routes
        for i, route_index in enumerate(tqdm(route_indices[:args.examples_per_cluster], desc=f"  Cluster {cluster_id} images")):
            if route_index >= len(all_routes):
                print(f"Warning: Route index {route_index} is out of bounds. Skipping.")
                continue
                
            # The route dict itself contains the full route data needed for drawing
            route_data = all_routes[route_index]
            
            image_path = cluster_dir / f"example_{i+1}_route_idx_{route_index}.{args.image_format}"
            generate_route_image(route_data, image_path)

    print("\n--- Visualization Complete ---")
    print(f"All cluster summaries and images have been saved to: {args.output_dir.resolve()}")
    print("------------------------------")

if __name__ == "__main__":
    main()