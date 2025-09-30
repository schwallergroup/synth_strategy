# clustering.py

import json
import os
import argparse
import numpy as np
import logging
import re
import ast
import glob
from pathlib import Path
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from typing import List, Dict, Any, Tuple

# Import the annotation function from our other module (package-relative)
from ..route_annotation.annotation import annotate_routes

# --- Standard Logging Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Core Clustering Helper Functions (largely unchanged, but with better docstrings) ---

def create_one_hot_vectors(routes: List[Dict[str, Any]]) -> Tuple[np.ndarray, List[str]]:
    """Creates one-hot encoded feature vectors from 'passing_functions' keys."""
    all_functions = set(func for route in routes for func in route.get('passing_functions', {}).keys())
    
    if not all_functions:
        logger.warning("No 'passing_functions' found. Cannot create vectors.")
        return np.array([]), []

    sorted_functions = sorted(list(all_functions))
    func_to_idx = {func: i for i, func in enumerate(sorted_functions)}

    vectors = np.zeros((len(routes), len(sorted_functions)), dtype=np.int8)
    for i, route in enumerate(routes):
        for func in route.get('passing_functions', {}).keys():
            if func in func_to_idx:
                vectors[i, func_to_idx[func]] = 1
    return vectors, sorted_functions

def find_optimal_clusters_silhouette(vectors: np.ndarray, max_k: int = 15) -> int:
    """Finds the optimal number of clusters (k) using the silhouette score."""
    k_range = range(2, min(vectors.shape[0], max_k + 1))
    if not list(k_range):
        logger.warning(f"Not enough samples ({vectors.shape[0]}) for silhouette analysis. Defaulting to 1 cluster.")
        return 1

    scores = []
    for k in k_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init='auto')
        labels = kmeans.fit_predict(vectors)
        scores.append(silhouette_score(vectors, labels))

    return k_range[np.argmax(scores)] if scores else 1

def load_function_metadata(code_dir: str) -> Dict[str, str]:
    """Load function descriptions from the metadata database JSON file.
    
    Args:
        code_dir: Directory where the metadata database might be located
    
    Returns:
        Dictionary mapping function filenames to their descriptions
    """
    metadata_map = {}
    
    # Try multiple possible locations for the metadata file
    # code_dir is typically "data/strategy_function_library/"
    code_path = Path(code_dir)
    possible_paths = [
        code_path / "function_metadata_database.json",  # In the same dir as functions
        code_path.parent / "function_metadata_database.json",  # In data/ directory
        code_path.parent / "3_metadata_output" / "function_metadata_database.json",
        code_path.parent.parent / "pipeline_output" / "3_metadata_output" / "function_metadata_database.json",
        Path("pipeline_output") / "3_metadata_output" / "function_metadata_database.json",
        Path("data") / "function_metadata_database.json",
    ]
    
    for metadata_path in possible_paths:
        if metadata_path.exists():
            try:
                with open(metadata_path, 'r') as f:
                    metadata_list = json.load(f)
                    for entry in metadata_list:
                        function_id = entry.get('function_id', '')
                        description = entry.get('description', '')
                        if function_id:
                            # Map both with and without .py extension
                            metadata_map[f"{function_id}.py"] = description
                            metadata_map[function_id] = description
                logger.info(f"Loaded metadata for {len(metadata_map)//2} functions from {metadata_path}")
                break
            except Exception as e:
                logger.warning(f"Could not load metadata from {metadata_path}: {e}")
    
    if not metadata_map:
        logger.warning("No function metadata database found. Descriptions will not be included.")
    
    return metadata_map

def get_cluster_defining_features(
    vectors: np.ndarray, 
    labels: np.ndarray, 
    model: KMeans, 
    names: List[str],
    metadata: Dict[str, str] = None
) -> Dict[str, List[Dict[str, Any]]]:
    """Identifies features that are most determinative for each cluster.
    
    Args:
        vectors: Feature vectors
        labels: Cluster labels
        model: Trained KMeans model
        names: Function names
        metadata: Optional dictionary mapping function names to descriptions
    
    Returns:
        Dictionary mapping cluster IDs to lists of defining features with scores and descriptions
    """
    defining_features = {}
    metadata = metadata or {}
    
    for i in range(model.n_clusters):
        in_mask = (labels == i)
        out_mask = ~in_mask
        p_in = model.cluster_centers_[i]
        p_out = np.mean(vectors[out_mask], axis=0) if np.any(out_mask) else np.zeros_like(p_in)
        distinctiveness = p_in - p_out
        
        # Create enriched feature list with descriptions
        features = []
        for name, score in sorted(zip(names, distinctiveness), key=lambda item: item[1], reverse=True):
            feature_entry = {
                "function": name,
                "score": float(score)
            }
            # Add description if available
            if name in metadata:
                feature_entry["description"] = metadata[name]
            features.append(feature_entry)
        
        defining_features[str(i)] = features
    
    return defining_features

def extract_docstrings(function_files: List[str], code_dir: str) -> Dict[str, str]:
    """Extracts docstrings for a list of function files from their source code."""
    docstrings = {}
    for func_file in function_files:
        file_path = Path(code_dir) / func_file
        if not file_path.exists():
            docstrings[func_file] = "Error: Source file not found."
            continue
        try:
            with open(file_path, "r") as f:
                tree = ast.parse(f.read())
            doc = "No docstring found."
            for node in ast.walk(tree):
                if isinstance(node, ast.FunctionDef) and node.name == "main":
                    doc = ast.get_docstring(node) or doc
                    break
            docstrings[func_file] = doc
        except Exception as e:
            docstrings[func_file] = f"Error parsing file: {e}"
    return docstrings

# --- Main API and Batch Entry Points ---

def perform_strategy_clustering(
    routes: List[Dict[str, Any]],
    code_dir: str,
    pre_annotated: bool = False,
    functions_dir: str = None
) -> Dict[str, Any]:
    """
    Performs strategy-based clustering on a list of synthesis routes.
    This function can either annotate routes first or use pre-annotated data.

    Args:
        routes: A list of synthesis route dictionaries.
        code_dir: Directory with source code files for extracting docstrings.
        pre_annotated: If True, assumes routes are already annotated. If False,
                       it will run annotation first (requires `functions_dir`).
        functions_dir: Directory with Python function files, required if `pre_annotated` is False.

    Returns:
        A dictionary containing the clustering results.
    """
    if not pre_annotated:
        if not functions_dir:
            raise ValueError("`functions_dir` is required when `pre_annotated` is False.")
        logger.info("Routes are not pre-annotated. Running annotation step first...")
        annotated_routes = annotate_routes(routes, functions_dir)
    else:
        logger.info("Using pre-annotated routes.")
        annotated_routes = routes

    # Use all routes for clustering (retrosynthesis tools output solved routes by default)
    logger.info(f"Processing {len(annotated_routes)} routes for clustering.")

    if len(annotated_routes) < 2:
        logger.warning("Not enough routes to perform clustering.")
        return {'status': 'skipped', 'reason': 'Not enough routes for clustering.'}

    # 1. Create feature vectors
    feature_vectors, unique_functions = create_one_hot_vectors(annotated_routes)
    if feature_vectors.shape[1] == 0:
        return {'status': 'skipped', 'reason': 'No passing functions found to create features.'}
    
    # 2. Extract docstrings for context
    function_docstrings = extract_docstrings(unique_functions, code_dir)
    
    # 3. Load function metadata descriptions
    function_metadata = load_function_metadata(code_dir)

    # 4. Find optimal K and perform final clustering
    optimal_k = find_optimal_clusters_silhouette(feature_vectors)
    kmeans = KMeans(n_clusters=optimal_k, random_state=42, n_init='auto')
    final_labels = kmeans.fit_predict(feature_vectors)

    # 5. Analyze and package results
    defining_features = get_cluster_defining_features(
        feature_vectors, final_labels, kmeans, unique_functions, function_metadata
    )
    
    # Map original route index to cluster id
    cluster_allocations = {str(i): int(label) for i, label in enumerate(final_labels)}

    return {
        'status': 'success',
        'optimal_k': optimal_k,
        'cluster_allocations': cluster_allocations,
        'cluster_defining_features': defining_features,
        'function_docstrings': function_docstrings
    }

def main():
    """
    Command-line interface for batch clustering. It aggregates all annotated
    JSON files from a directory and performs clustering on the combined dataset.
    """
    parser = argparse.ArgumentParser(description="Aggregate annotated files and perform strategy clustering.")
    parser.add_argument("--input-dir", required=True, help="Directory containing ANNOTATED JSON files (outputs from annotation.py).")
    parser.add_argument("--output", required=True, help="Path to save the final clustering analysis JSON.")
    parser.add_argument("--code-dir", required=True, help="Path to the source code files for functions to extract docstrings.")
    args = parser.parse_args()

    # Aggregate all routes from annotated files
    all_annotated_routes = []
    json_files = glob.glob(os.path.join(args.input_dir, '*.json'))
    logger.info(f"Found {len(json_files)} annotated files to process in {args.input_dir}.")
    for file_path in json_files:
        try:
            with open(file_path, 'r') as f:
                all_annotated_routes.extend(json.load(f))
        except (json.JSONDecodeError, IOError) as e:
            logger.warning(f"Could not load or parse {file_path}: {e}")

    # Perform clustering on the aggregated, pre-annotated data
    results = perform_strategy_clustering(
        routes=all_annotated_routes,
        code_dir=args.code_dir,
        pre_annotated=True  # Data is already annotated
    )

    # Save final analysis
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2)
    logger.info(f"Clustering analysis complete. Results saved to {args.output}")

if __name__ == "__main__":
    main()