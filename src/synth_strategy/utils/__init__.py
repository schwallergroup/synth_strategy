"""
Utility functions for SynthStrategy.

This module provides common utility functions used across the SynthStrategy package.
"""

import logging
import sys
from pathlib import Path
from typing import List, Dict, Any, Optional
import json


def setup_logging(verbose: bool = False, log_level: Optional[str] = None) -> None:
    """
    Setup logging configuration.
    
    Args:
        verbose: Enable verbose logging
        log_level: Specific log level to use
    """
    if log_level is None:
        log_level = "DEBUG" if verbose else "INFO"
    
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout)
        ]
    )


def validate_paths(paths: Dict[str, str], must_exist: bool = True) -> None:
    """
    Validate that paths exist and are accessible.
    
    Args:
        paths: Dictionary of path names to path strings
        must_exist: Whether paths must exist
        
    Raises:
        FileNotFoundError: If required paths don't exist
        ValueError: If paths are invalid
    """
    for name, path_str in paths.items():
        if not path_str:
            raise ValueError(f"{name} path is empty")
        
        path = Path(path_str)
        
        if must_exist and not path.exists():
            raise FileNotFoundError(f"{name} path does not exist: {path}")
        
        if must_exist and not path.is_file() and not path.is_dir():
            raise ValueError(f"{name} path is neither file nor directory: {path}")


def load_json_safe(file_path: str, default: Any = None) -> Any:
    """
    Safely load JSON file with error handling.
    
    Args:
        file_path: Path to JSON file
        default: Default value to return if loading fails
        
    Returns:
        Loaded JSON data or default value
    """
    try:
        with open(file_path, 'r') as f:
            return json.load(f)
    except (json.JSONDecodeError, IOError) as e:
        logging.warning(f"Could not load JSON file {file_path}: {e}")
        return default


def save_json_safe(data: Any, file_path: str) -> bool:
    """
    Safely save data to JSON file with error handling.
    
    Args:
        data: Data to save
        file_path: Path to save file
        
    Returns:
        True if successful, False otherwise
    """
    try:
        Path(file_path).parent.mkdir(parents=True, exist_ok=True)
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2)
        return True
    except (IOError, TypeError) as e:
        logging.error(f"Could not save JSON file {file_path}: {e}")
        return False


def find_json_files(directory: str, pattern: str = "*.json") -> List[str]:
    """
    Find all JSON files in a directory.
    
    Args:
        directory: Directory to search
        pattern: File pattern to match
        
    Returns:
        List of file paths
    """
    import glob
    
    search_path = str(Path(directory) / pattern)
    return glob.glob(search_path)


def ensure_directory(path: str) -> Path:
    """
    Ensure directory exists, creating it if necessary.
    
    Args:
        path: Directory path
        
    Returns:
        Path object
    """
    dir_path = Path(path)
    dir_path.mkdir(parents=True, exist_ok=True)
    return dir_path


def get_file_size_mb(file_path: str) -> float:
    """
    Get file size in megabytes.
    
    Args:
        file_path: Path to file
        
    Returns:
        File size in MB
    """
    return Path(file_path).stat().st_size / (1024 * 1024)


def format_duration(seconds: float) -> str:
    """
    Format duration in seconds to human-readable string.
    
    Args:
        seconds: Duration in seconds
        
    Returns:
        Formatted duration string
    """
    if seconds < 60:
        return f"{seconds:.1f}s"
    elif seconds < 3600:
        minutes = seconds / 60
        return f"{minutes:.1f}m"
    else:
        hours = seconds / 3600
        return f"{hours:.1f}h"


def progress_bar(iterable, desc: str = "Processing", total: Optional[int] = None):
    """
    Create a progress bar for iterables.
    
    Args:
        iterable: Iterable to wrap
        desc: Description for progress bar
        total: Total number of items (if known)
        
    Yields:
        Items from iterable
    """
    try:
        from tqdm import tqdm
        return tqdm(iterable, desc=desc, total=total)
    except ImportError:
        # Fallback if tqdm is not available
        return iterable


def check_dependencies() -> Dict[str, bool]:
    """
    Check if required dependencies are available.
    
    Returns:
        Dictionary mapping dependency names to availability
    """
    dependencies = {
        "numpy": False,
        "scikit-learn": False,
        "rdkit": False,
        "aizynthfinder": False,
        "tqdm": False,
        "yaml": False,
    }
    
    for dep in dependencies:
        try:
            if dep == "yaml":
                import yaml
            elif dep == "scikit-learn":
                import sklearn
            else:
                __import__(dep)
            dependencies[dep] = True
        except ImportError:
            dependencies[dep] = False
    
    return dependencies


def print_dependency_status() -> None:
    """Print status of required dependencies."""
    deps = check_dependencies()
    
    print("Dependency Status:")
    print("-" * 20)
    
    for dep, available in deps.items():
        status = "✅" if available else "❌"
        print(f"{status} {dep}")
    
    missing = [dep for dep, available in deps.items() if not available]
    if missing:
        print(f"\nMissing dependencies: {', '.join(missing)}")
        print("Install with: pip install " + " ".join(missing))


__all__ = [
    'setup_logging',
    'validate_paths',
    'load_json_safe',
    'save_json_safe',
    'find_json_files',
    'ensure_directory',
    'get_file_size_mb',
    'format_duration',
    'progress_bar',
    'check_dependencies',
    'print_dependency_status',
]