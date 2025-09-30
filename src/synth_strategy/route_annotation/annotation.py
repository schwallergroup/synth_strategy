# annotation.py

import sys
import json
import argparse
import logging
import importlib.util
import inspect
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Any
import os

# --- Standard Logging Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

# --- Global variables for worker processes ---
WORKER_FUNCTIONS = {}

# --- Helper Functions for Route Pre-processing ---

def traverse_and_flip(node: dict):
    """Recursively flips 'mapped_reaction_smiles' from retro to forward direction."""
    if node.get("type") == "reaction":
        metadata = node.get("metadata", {})
        mapped_smiles = metadata.get("mapped_reaction_smiles")
        if isinstance(mapped_smiles, str) and ">>" in mapped_smiles:
            target, precursors = mapped_smiles.split(">>", 1)
            node["metadata"]["mapped_reaction_smiles"] = f"{precursors}>>{target}"

    if "children" in node and isinstance(node["children"], list):
        for child_node in node["children"]:
            traverse_and_flip(child_node)

def calculate_depth(node: Dict[str, Any]) -> int:
    """Calculate the maximum depth of a synthesis route tree."""
    if not node.get('children'):
        return 0
    return 1 + max(calculate_depth(child) for child in node['children']) if node.get('children') else 0

def sort_route_by_depth(route: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively sort each node's children so that the deepest subtree is first."""
    sorted_route = route.copy()
    if 'children' in sorted_route and sorted_route['children']:
        sorted_children = [sort_route_by_depth(child) for child in sorted_route['children']]
        sorted_children.sort(key=calculate_depth, reverse=True)
        sorted_route['children'] = sorted_children
    return sorted_route


# --- Worker Initialization and Logic ---

def init_worker(directory_path: str, function_name: str):
    """Initializer for each worker process. Loads functions from files once."""
    global WORKER_FUNCTIONS
    logger.info(f"Worker process (PID: {os.getpid()}) initializing... Loading functions from {directory_path}")
    directory = Path(directory_path)
    if not directory.is_dir():
        return

    for py_file in directory.glob("*.py"):
        try:
            spec = importlib.util.spec_from_file_location(py_file.stem, py_file)
            if spec and spec.loader:
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                func = getattr(module, function_name, None)
                if inspect.isfunction(func):
                    WORKER_FUNCTIONS[py_file.name] = func
        except Exception as e:
            logger.error(f"[Worker] Failed to load {py_file.name}: {e}")
    logger.info(f"[Worker] Initialization complete. Loaded {len(WORKER_FUNCTIONS)} functions.")

def process_route_in_worker(route: Dict[str, Any]) -> Dict[str, Any]:
    """Applies all loaded functions to a single route and annotates it with results."""
    # Pre-process the route: sort by depth and flip reaction SMILES
    processed_route = sort_route_by_depth(route)
    #traverse_and_flip(processed_route)

    passing_functions = {}  # Store results as {filename: findings}
    errored_functions = []

    for filename, func in WORKER_FUNCTIONS.items():
        try:
            result = func(processed_route.copy()) # Pass copy to prevent side-effects

            # Handle both (bool, findings) and simple bool return types
            if isinstance(result, tuple) and len(result) == 2:
                is_passing, findings = result
                if is_passing:
                    passing_functions[filename] = findings
            elif isinstance(result, bool) and result:
                passing_functions[filename] = [] # Default to empty findings
            
        except Exception as e:
            logger.error(f"Error in {filename}: {e}")
            errored_functions.append(filename)

    processed_route['passing_functions'] = passing_functions
    processed_route['errored_functions'] = sorted(errored_functions)
    return processed_route

# --- Main Entry Points (API and Batch) ---

def annotate_routes(routes: List[Dict[str, Any]], functions_dir: str, num_workers: int = None) -> List[Dict[str, Any]]:
    """
    Annotates a list of synthesis routes by applying a set of functions in parallel.
    This is the primary function to be called by an API.

    Args:
        routes: A list of synthesis route dictionaries.
        functions_dir: The directory containing the Python function files.
        num_workers: The number of processes to use. Defaults to the number of CPU cores.

    Returns:
        A list of the same routes, now annotated with 'passing_functions' and 'errored_functions'.
    """
    if not routes:
        return []
        
    if num_workers is None:
        num_workers = os.cpu_count()

    logger.info(f"Annotating {len(routes)} routes using up to {num_workers} workers...")
    annotated_routes = []

    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=init_worker,
        initargs=(functions_dir, 'main')
    ) as executor:
        
        future_to_route = {executor.submit(process_route_in_worker, route): route for route in routes}
        
        for i, future in enumerate(as_completed(future_to_route), 1):
            try:
                result = future.result()
                annotated_routes.append(result)
                if i % 10 == 0 or i == len(routes):
                     logger.info(f"({i}/{len(routes)}) routes annotated.")
            except Exception as e:
                logger.error(f"A critical error occurred while processing a route: {e}", exc_info=True)
    
    return annotated_routes

def main():
    """
    Command-line interface for batch processing a single JSON file of routes.
    Designed to be called by a SLURM job.
    """
    parser = argparse.ArgumentParser(description='Apply strategy functions to a JSON file of synthesis routes.')
    parser.add_argument('--input', type=str, required=True, help='Path to the input JSON file containing routes.')
    parser.add_argument('--output', type=str, required=True, help='Path to save the annotated output JSON file.')
    parser.add_argument('--functions-dir', type=str, required=True, help='Directory containing Python function files.')
    parser.add_argument('--workers', type=int, default=4, help='Number of worker processes.')
    args = parser.parse_args()

    # Load data from input file
    try:
        with open(args.input, 'r') as f:
            raw_data = json.load(f)

        # Handle both list-of-routes and the HDF5-like format
        if isinstance(raw_data, dict) and 'data' in raw_data and raw_data.get('data'):
            trees_data = raw_data['data'][0].get('trees')
            all_routes = json.loads(trees_data) if isinstance(trees_data, str) else trees_data
        elif isinstance(raw_data, list):
            all_routes = raw_data
        else:
            raise ValueError("Unsupported JSON structure.")
    except Exception as e:
        logger.error(f"Failed to load and parse dataset from {args.input}: {e}")
        sys.exit(1)

    # Annotate the routes
    final_routes = annotate_routes(
        routes=all_routes,
        functions_dir=args.functions_dir,
        num_workers=args.workers
    )

    # Save the results
    logger.info(f"Saving {len(final_routes)} annotated routes to {args.output}...")
    try:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(args.output, 'w') as f:
            json.dump(final_routes, f, indent=2)
        logger.info("Annotation complete.")
    except Exception as e:
        logger.error(f"Failed to save results to {args.output}: {e}")

if __name__ == "__main__":
    main()