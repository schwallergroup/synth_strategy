#!/usr/bin/env python3
"""
Test script for strategy annotation and clustering using the SynthStrategyAPI.

This script loads the first 5 routes from data/routes/val_t.json and tests:
1. Strategy annotation functionality
2. Strategy clustering functionality

Usage:
    python test.py
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Any

# Determine the correct path structure
current_dir = Path(__file__).parent
possible_paths = [
    current_dir / "src" / "synth_strategy",  # If running from project root
    current_dir / "synth_strategy",          # If src is current directory
    current_dir,                             # If synth_strategy is current directory
]

synth_strategy_path = None
for path in possible_paths:
    if (path / "api.py").exists():
        synth_strategy_path = path
        break

if synth_strategy_path is None:
    print("Could not find synth_strategy directory with api.py")
    print("Searched in:")
    for path in possible_paths:
        print(f"  {path}")
    sys.exit(1)

print(f"Found synth_strategy at: {synth_strategy_path}")

# --- FIX 1: Add the PARENT directory ('src') to the path, not the package directory itself ---
project_src_path = synth_strategy_path.parent
sys.path.insert(0, str(project_src_path))
print(f"Added to sys.path: {project_src_path}")

try:
    # --- FIX 2: Use the full, absolute import from the package name ---
    from synth_strategy.api import SynthStrategyAPI
except ImportError as e:
    print(f"Error importing SynthStrategyAPI: {e}")
    # The error message from before might show up here if something is still wrong.
    # The new path added is the parent of your package.
    print(f"Tried to import from: {project_src_path}") 
    print("Make sure the synth_strategy directory structure is correct.")
    print("Available items in the added path:")
    if project_src_path.exists():
        for item in project_src_path.iterdir():
            print(f"  {item.name}")
    sys.exit(1)

# (The rest of your test.py file remains the same)

def load_test_data(data_path: str = "data/routes/val_t.json", num_routes: int = 5) -> List[Dict[str, Any]]:
    """
    Load the first N routes from the validation dataset.
    
    Args:
        data_path: Path to the validation data JSON file
        num_routes: Number of routes to load for testing
        
    Returns:
        List of route dictionaries
    """
    try:
        with open(data_path, 'r') as f:
            all_routes = json.load(f)
        
        if not isinstance(all_routes, list):
            raise ValueError(f"Expected a list in {data_path}, got {type(all_routes)}")
        
        if len(all_routes) < num_routes:
            print(f"Warning: Only {len(all_routes)} routes available, using all of them.")
            num_routes = len(all_routes)
        
        test_routes = all_routes[:num_routes]
        print(f"Loaded {len(test_routes)} routes from {data_path}")
        return test_routes
        
    except FileNotFoundError:
        print(f"Error: Could not find {data_path}")
        print("Make sure you're running this script from the project root directory.")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error parsing JSON from {data_path}: {e}")
        sys.exit(1)


def print_route_summary(routes: List[Dict[str, Any]]) -> None:
    """Print a summary of the loaded routes."""
    print("\n" + "="*50)
    print("ROUTE SUMMARY")
    print("="*50)
    
    for i, route in enumerate(routes):
        print(f"\nRoute {i+1}:")
        # Print some basic info about each route
        route_keys = list(route.keys())
        print(f"  Keys: {route_keys}")
        
        # Try to show some identifying information
        if 'id' in route:
            print(f"  ID: {route['id']}")
        if 'target' in route:
            print(f"  Target: {str(route['target'])[:100]}...")
        if 'route' in route:
            print(f"  Route length: {len(route['route']) if isinstance(route['route'], list) else 'N/A'}")


def test_annotation(api: SynthStrategyAPI, routes: List[Dict[str, Any]], functions_dir: str) -> List[Dict[str, Any]]:
    """
    Test the strategy annotation functionality.
    
    Args:
        api: The SynthStrategyAPI instance
        routes: List of routes to annotate
        functions_dir: Path to strategy functions directory
        
    Returns:
        Annotated routes
    """
    print("\n" + "="*50)
    print("TESTING STRATEGY ANNOTATION")
    print("="*50)
    
    # Check if functions directory exists
    functions_path = Path(functions_dir)
    if not functions_path.exists():
        print(f"Error: Strategy functions directory '{functions_dir}' not found.")
        print("Please specify the correct path to your strategy functions.")
        return routes
    
    # Count Python files in the directory
    py_files = list(functions_path.glob("*.py"))
    print(f"Found {len(py_files)} Python files in {functions_dir}")
    
    try:
        # Run annotation
        print(f"Annotating {len(routes)} routes...")
        annotated_routes = api.annotate_strategies(
            routes=routes,
            functions_dir=functions_dir,
            num_workers=2  # Use fewer workers for testing
        )
        
        print("\nAnnotation Results:")
        print("-" * 30)
        
        for i, route in enumerate(annotated_routes):
            print(f"\nRoute {i+1}:")
            
            # Check for annotation results
            passing_funcs = route.get('passing_functions', [])
            errored_funcs = route.get('errored_functions', [])
            
            print(f"  Passing functions: {len(passing_funcs)}")
            if passing_funcs:
                print(f"    {passing_funcs[:3]}{'...' if len(passing_funcs) > 3 else ''}")
            
            print(f"  Errored functions: {len(errored_funcs)}")
            if errored_funcs:
                print(f"    {list(errored_funcs.keys())[:3]}{'...' if len(errored_funcs) > 3 else ''}")
        
        return annotated_routes
        
    except Exception as e:
        print(f"Error during annotation: {e}")
        print("Returning original routes without annotation.")
        return routes


def test_clustering(api: SynthStrategyAPI, annotated_routes: List[Dict[str, Any]], code_dir: str) -> Dict[str, Any]:
    """
    Test the strategy clustering functionality.
    
    Args:
        api: The SynthStrategyAPI instance
        annotated_routes: List of annotated routes
        code_dir: Path to strategy functions source code directory
        
    Returns:
        Clustering results dictionary
    """
    print("\n" + "="*50)
    print("TESTING STRATEGY CLUSTERING")
    print("="*50)
    
    # Check if any routes have annotations
    has_annotations = any(
        route.get('passing_functions') or route.get('errored_functions') 
        for route in annotated_routes
    )
    
    if not has_annotations:
        print("Warning: No annotated routes found. Clustering may not work properly.")
        print("Make sure the annotation step completed successfully.")
    
    try:
        # Run clustering
        print(f"Clustering {len(annotated_routes)} annotated routes...")
        clustering_results = api.cluster_strategies(
            annotated_routes=annotated_routes,
            code_dir=code_dir
        )
        
        print("\nClustering Results:")
        print("-" * 30)
        
        # Print key results
        if 'optimal_k' in clustering_results:
            print(f"Optimal number of clusters: {clustering_results['optimal_k']}")
        
        if 'cluster_assignments' in clustering_results:
            assignments = clustering_results['cluster_assignments']
            print(f"Cluster assignments: {assignments}")
        
        if 'cluster_summaries' in clustering_results:
            summaries = clustering_results['cluster_summaries']
            print(f"Found {len(summaries)} cluster summaries")
            for i, summary in enumerate(summaries[:3]):  # Show first 3
                print(f"  Cluster {i}: {str(summary)[:100]}...")
        
        # Print any other interesting keys
        other_keys = [k for k in clustering_results.keys() 
                     if k not in ['optimal_k', 'cluster_assignments', 'cluster_summaries']]
        if other_keys:
            print(f"Other result keys: {other_keys}")
        
        return clustering_results
        
    except Exception as e:
        print(f"Error during clustering: {e}")
        return {}


def main():
    """Main test function."""
    print("Starting annotation and clustering test...")
    
    # Configuration - adjust these paths as needed
    # Try to find the data file in various locations
    possible_data_paths = [
        "data/routes/val_t.json",
        "../data/routes/val_t.json",
        "../../data/routes/val_t.json",
    ]
    
    DATA_PATH = None
    for path in possible_data_paths:
        if Path(path).exists():
            DATA_PATH = path
            break
    
    if DATA_PATH is None:
        print("Could not find val_t.json in any of the expected locations:")
        for path in possible_data_paths:
            print(f"  {path}")
        print("Please adjust the DATA_PATH in the script or create some test data.")
        return
    
    # Look for strategy functions directory
    possible_func_dirs = [
        "strategy_functions",
        "src/synth_strategy/strategy_functions", 
        "synth_strategy/strategy_functions",
        "functions",
    ]
    
    FUNCTIONS_DIR = None
    for func_dir in possible_func_dirs:
        if Path(func_dir).exists():
            FUNCTIONS_DIR = func_dir
            break
    
    if FUNCTIONS_DIR is None:
        print("Could not find strategy functions directory. Looked in:")
        for path in possible_func_dirs:
            print(f"  {path}")
        print("Creating a dummy functions directory for testing...")
        FUNCTIONS_DIR = "test_functions"
        Path(FUNCTIONS_DIR).mkdir(exist_ok=True)
        # Create a simple test function
        test_func = Path(FUNCTIONS_DIR) / "test_strategy.py"
        test_func.write_text('''
def test_strategy(route_dict):
    """A simple test strategy function."""
    return True
''')
    
    CODE_DIR = FUNCTIONS_DIR
    NUM_TEST_ROUTES = 5
    
    # Load test data
    routes = load_test_data(DATA_PATH, NUM_TEST_ROUTES)
    print_route_summary(routes)
    
    # Initialize API
    try:
        api = SynthStrategyAPI()
    except Exception as e:
        print(f"Error initializing SynthStrategyAPI: {e}")
        sys.exit(1)
    
    # Test annotation
    annotated_routes = test_annotation(api, routes, FUNCTIONS_DIR)
    
    # Test clustering
    clustering_results = test_clustering(api, annotated_routes, CODE_DIR)
    
    print("\n" + "="*50)
    print("TEST SUMMARY")
    print("="*50)
    print(f"✓ Loaded {len(routes)} test routes")
    print(f"✓ Annotation {'completed' if any('passing_functions' in r for r in annotated_routes) else 'attempted'}")
    print(f"✓ Clustering {'completed' if clustering_results else 'attempted'}")
    
    if clustering_results and 'optimal_k' in clustering_results:
        print(f"✓ Found {clustering_results['optimal_k']} optimal clusters")
    
    print("\nTest completed!")


if __name__ == "__main__":
    main()