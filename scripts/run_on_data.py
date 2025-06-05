import sys
import json
import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict
from steerable_retro.api import StrategyDescriptor

def process_route(route: Dict) -> Dict:
    """Process a single route using StrategyDescriptor."""
    descriptor = StrategyDescriptor()
    result = descriptor(route)
    descriptor.results_stats(result)
    return result

def process_dataset(file_path: str, num_workers: int = -1) -> List[Dict]:
    """
    Process the dataset with parallel execution.
    
    Args:
        file_path: Path to the JSON file containing routes
        num_workers: Number of worker processes (-1 for CPU count)
    
    Returns:
        List of processed route descriptors
    """
    # Read the JSON file
    try:
        with open(file_path, 'r') as f:
            routes = json.load(f)
    except json.JSONDecodeError:
        print(f"Error: Invalid JSON file: {file_path}")
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}")
        sys.exit(1)

    # Ensure routes is a list
    if not isinstance(routes, list):
        routes = [routes]

    results = []
    
    # Process routes in parallel
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        # Submit all routes for processing
        future_to_route = {executor.submit(process_route, route): route for route in routes}
        
        # Collect results as they complete
        for future in as_completed(future_to_route):
            try:
                result = future.result()
                results.append(result)
            except Exception as e:
                print(f"Error processing route: {e}")

    return results

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Process routes using StrategyDescriptor')
    parser.add_argument('dataset', type=str, help='Path to the JSON dataset file')
    parser.add_argument('--workers', type=int, default=-1,
                       help='Number of worker processes (-1 for CPU count)')
    parser.add_argument('--output', type=str, help='Output file path for results (optional)')
    return parser.parse_args()

def main():
    """Main function."""
    # Parse arguments
    args = parse_arguments()
    
    # Process the dataset
    print(f"Processing dataset: {args.dataset}")
    results = process_dataset(args.dataset, args.workers)
    print(f"Processed {len(results)} routes")

    # Save results if output path is specified
    if args.output:
        output_path = Path(args.output)
        try:
            with open(output_path, 'w') as f:
                json.dump(results, f, indent=2)
            print(f"Results saved to: {output_path}")
        except Exception as e:
            print(f"Error saving results: {e}")

    return results

if __name__ == "__main__":
    main()
