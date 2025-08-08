import importlib.util
import inspect
import logging
from pathlib import Path
from typing import Dict, List, Callable, Tuple

# It's good practice to have a logger
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class StrategyDescriptor:
    """
    Loads all functions from a directory once and then applies them to routes.
    """
    def __init__(self, directory_path: str, function_name: str = "main"):
        """
        Loads all functions ONCE upon initialization.

        Args:
            directory_path: Directory containing the Python strategy files.
            function_name: The name of the function to execute in each file (e.g., "main").
        """
        self.functions: Dict[str, Callable] = {}
        self.load_errors: List[str] = []
        self._load_functions_from_directory(directory_path, function_name)

    def _load_functions_from_directory(self, directory_path: str, function_name: str):
        """Private method to load functions and track their source files."""
        directory = Path(directory_path)
        if not directory.is_dir():
            logger.error(f"Directory does not exist: {directory_path}")
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
                        self.functions[py_file.name] = func
                    else:
                        logger.warning(f"Could not find function '{function_name}' in {py_file.name}")
                        self.load_errors.append(py_file.name)
                
            except Exception as e:
                # This catches syntax errors or other import-time problems
                logger.error(f"Failed to load {py_file.name} due to error: {e}")
                self.load_errors.append(py_file.name)
        
        logger.info(f"Successfully loaded {len(self.functions)} functions.")
        if self.load_errors:
            logger.warning(f"Failed to load or find function in {len(self.load_errors)} files: {self.load_errors}")

    def __call__(self, route: Dict) -> Tuple[List[str], List[str]]:
        """
        Executes all loaded functions on a single route.

        Args:
            route: A single route dictionary.

        Returns:
            A tuple containing two lists:
            - passing_files: List of filenames whose function returned True.
            - errored_files: List of filenames whose function raised an exception.
        """
        passing_files = []
        errored_files = []

        for filename, func in self.functions.items():
            try:
                # Run the function on the route
                if func(route.copy()) is True:
                    passing_files.append(filename)
            except Exception as e:
                # If the function itself errors during execution
                logger.debug(f"Function from '{filename}' raised an error on route: {e}", exc_info=False)
                errored_files.append(filename)
        
        return passing_files, errored_files

def add_metadata_to_route(route: Dict, passing_files: List[str], errored_files: List[str]) -> Dict:
    """
    Adds 'passing_functions' and 'errored_functions' keys to the top-level node of a route.
    """
    modified_route = route.copy()
    modified_route['passing_functions'] = passing_files
    modified_route['errored_functions'] = errored_files
    return modified_route