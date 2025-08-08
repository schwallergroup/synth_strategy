import sys
import json
import argparse
import logging
import importlib.util
import inspect
import traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Any, Tuple
from io import StringIO
from contextlib import redirect_stdout, redirect_stderr

# --- Global variables for worker processes ---
# These will be populated by the init_worker function once in each worker.
WORKER_FUNCTIONS = {}
WORKER_LOAD_ERRORS = []

# --- Standard Logging Setup ---
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def init_worker(directory_path: str, function_name: str):
    """
    Initializer function for each worker process. Loads all functions from files once.
    This function is run by the ProcessPoolExecutor for each new worker process.
    """
    global WORKER_FUNCTIONS, WORKER_LOAD_ERRORS
    
    logger.info(f"Worker process initializing... Loading functions from {directory_path}")
    
    directory = Path(directory_path)
    if not directory.is_dir():
        logger.error(f"[Worker] Directory does not exist: {directory_path}")
        return

    for py_file in directory.glob("*.py"):
        try:
            spec = importlib.util.spec_from_file_location(py_file.stem, py_file)
            if spec and spec.loader:
                module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(module)
                
                func = getattr(module, function_name, None)
                if inspect.isfunction(func):
                    # Store the function with its filename as the key
                    WORKER_FUNCTIONS[py_file.name] = func
                else:
                    logger.warning(f"[Worker] Could not find function '{function_name}' in {py_file.name}")
                    WORKER_LOAD_ERRORS.append(py_file.name)
            
        except Exception as e:
            # This catches syntax errors or other import-time problems
            logger.error(f"[Worker] Failed to load {py_file.name} due to error: {e}")
            WORKER_LOAD_ERRORS.append(py_file.name)
    
    logger.info(f"[Worker] Initialization complete. Loaded {len(WORKER_FUNCTIONS)} functions.")

def process_route_in_worker(route: Dict[str, Any]) -> Dict[str, Any]:
    """
    This is the target function for each worker. It runs on a single route
    and uses the pre-loaded global functions.
    """
    passing_files = []
    errored_files = []
    function_debug_info = {}

    for filename, func in WORKER_FUNCTIONS.items():
        # Capture stdout and stderr
        stdout_capture = StringIO()
        stderr_capture = StringIO()
        
        debug_entry = {
            'stdout': '',
            'stderr': '',
            'error': None,
            'traceback': None,
            'result': None,
            'status': 'unknown'
        }
        
        try:
            # IMPORTANT: The provided code_XX.py files modify data structures in place.
            # Passing a copy prevents one function from affecting the input for the next.
            with redirect_stdout(stdout_capture), redirect_stderr(stderr_capture):
                result = func(route.copy())
            
            # Capture the output
            debug_entry['stdout'] = stdout_capture.getvalue()
            debug_entry['stderr'] = stderr_capture.getvalue()
            debug_entry['result'] = result
            
            if result is True:
                passing_files.append(filename)
                debug_entry['status'] = 'passed'
            else:
                debug_entry['status'] = 'failed'
                
        except Exception as e:
            # If the function itself errors during execution on the route
            errored_files.append(filename)
            debug_entry['error'] = str(e)
            debug_entry['traceback'] = traceback.format_exc()
            debug_entry['status'] = 'error'
            debug_entry['stdout'] = stdout_capture.getvalue()
            debug_entry['stderr'] = stderr_capture.getvalue()
        
        function_debug_info[filename] = debug_entry

    # Add the metadata to the top-level node of the route
    route['passing_functions'] = sorted(passing_files)
    route['errored_functions'] = sorted(errored_files)
    route['function_debug_info'] = function_debug_info
    # You can also store files that failed to load if you want
    route['load_errors'] = WORKER_LOAD_ERRORS.copy() if WORKER_LOAD_ERRORS else []
    
    return route

def save_debug_info(processed_routes: List[Dict[str, Any]], debug_output_path: str):
    """
    Extract and save debug information to a separate file for easier analysis.
    """
    debug_data = []
    
    for idx, route in enumerate(processed_routes):
        route_id = route.get('smiles', route.get('id', 'unknown')) + f"_{idx}"
        
        route_debug = {
            'route_id': route_id,
            'load_errors': route.get('load_errors', []),
            'function_results': {}
        }
        
        function_debug_info = route.get('function_debug_info', {})
        for filename, debug_info in function_debug_info.items():
            route_debug['function_results'][filename] = debug_info
        
        debug_data.append(route_debug)
    
    try:
        with open(debug_output_path, 'w') as f:
            json.dump(debug_data, f, indent=2)
        logger.info(f"Debug information saved to {debug_output_path}")
    except Exception as e:
        logger.error(f"Failed to save debug information to {debug_output_path}: {e}")

def run_parallel_processing(file_path: str, functions_directory: str, function_name: str, num_workers: int, limit: int):
    """
    Main orchestrator for processing the dataset in parallel.
    """
    try:
        with open(file_path, 'r') as f:
            all_routes = json.load(f)
    except (json.JSONDecodeError, FileNotFoundError) as e:
        logger.error(f"Failed to load dataset from {file_path}: {e}")
        sys.exit(1)

    routes_to_process = all_routes if limit == -1 else all_routes[:limit]
    logger.info(f"Preparing to process {len(routes_to_process)} routes using {num_workers} workers...")
    
    processed_routes = []
    
    # Use ProcessPoolExecutor with an initializer
    with ProcessPoolExecutor(
        max_workers=num_workers,
        initializer=init_worker,
        initargs=(functions_directory, function_name)
    ) as executor:
        
        # Submit jobs
        future_to_route = {
            executor.submit(process_route_in_worker, route): route
            for route in routes_to_process
        }
        
        for i, future in enumerate(as_completed(future_to_route), 1):
            original_route = future_to_route[future]
            try:
                # The result is the route with 'passing_functions' and 'errored_functions' added
                annotated_route = future.result()
                processed_routes.append(annotated_route)
                logger.info(f"({i}/{len(routes_to_process)}) Processed route: {original_route.get('smiles', 'N/A')}")
            except Exception as e:
                # This catches errors in the worker logic itself, not the strategy functions
                logger.error(f"A critical error occurred while processing route {original_route.get('smiles', 'N/A')}: {e}", exc_info=True)

    return processed_routes

def main():
    parser = argparse.ArgumentParser(description='Process chemical routes and annotate them with strategy function results.')
    parser.add_argument('dataset', type=str, help='Path to the JSON dataset file containing routes.')
    parser.add_argument('--functions-dir', type=str, required=True, help='Directory containing Python function files.')
    parser.add_argument('--function-name', type=str, default='main', help='Name of the function to execute in each file (default: main).')
    parser.add_argument('--workers', type=int, default=20, help='Number of worker processes.')
    parser.add_argument('--output', type=str, required=True, help='Output file path for the annotated routes.')
    parser.add_argument('--debug-output', type=str, help='Output file path for debug information (default: <output>_debug.json).')
    parser.add_argument('--limit', type=int, default=1000000, help='Limit the number of routes to process for testing (default: all).')
    args = parser.parse_args()

    # Set default debug output path if not provided
    if not args.debug_output:
        output_path = Path(args.output)
        args.debug_output = str(output_path.parent / f"{output_path.stem}_debug.json")

    # Run the processing
    final_routes = run_parallel_processing(
        file_path=args.dataset,
        functions_directory=args.functions_dir,
        function_name=args.function_name,
        num_workers=args.workers,
        limit=args.limit
    )

    # Save the debug information
    logger.info(f"Saving debug information to {args.debug_output}...")
    save_debug_info(final_routes, args.debug_output)

    # Clean the routes for the main output (remove debug info to keep it clean)
    for route in final_routes:
        route.pop('function_debug_info', None)
        route.pop('load_errors', None)

    # Save the main results
    logger.info(f"Saving {len(final_routes)} processed routes to {args.output}...")
    try:
        with open(args.output, 'w') as f:
            json.dump(final_routes, f, indent=2)
        logger.info("Processing complete.")
    except Exception as e:
        logger.error(f"Failed to save results to {args.output}: {e}")

if __name__ == "__main__":
    main()