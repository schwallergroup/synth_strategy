"""
Quick access functions for SynthStrategy.

This module provides simplified, one-liner functions for common operations
including clustering, retrieval, and basic analysis.
"""

import asyncio
from pathlib import Path
from typing import List, Dict, Any, Optional, Union

from .api import SynthStrategyAPI
from .config import load_config, Config
from .utils import load_json_safe, save_json_safe, find_json_files


def quick_cluster(
    routes_path: Union[str, List[Dict[str, Any]]],
    functions_dir: str,
    output_dir: Optional[str] = None,
    config_path: Optional[str] = None,
    visualize: bool = True
) -> Dict[str, Any]:
    """
    Quick clustering of synthesis routes.
    
    Args:
        routes_path: Path to routes directory/file or list of route dictionaries
        functions_dir: Directory containing strategy functions
        output_dir: Output directory for results
        config_path: Path to configuration file
        visualize: Whether to generate visualizations
        
    Returns:
        Clustering results dictionary
    """
    # Load configuration
    config = load_config(config_path) if config_path else Config()
    
    # Initialize API
    api = SynthStrategyAPI()
    
    # Load routes
    if isinstance(routes_path, str):
        if Path(routes_path).is_dir():
            routes = _load_routes_from_directory(routes_path)
        else:
            routes = load_json_safe(routes_path, [])
    else:
        routes = routes_path
    
    if not routes:
        raise ValueError("No routes found to cluster")
    
    print(f"🔬 Clustering {len(routes)} routes...")
    
    # Perform clustering
    results = api.cluster_strategies(
        annotated_routes=routes,
        code_dir=functions_dir
    )
    
    if results.get('status') != 'success':
        raise RuntimeError(f"Clustering failed: {results.get('reason', 'Unknown error')}")
    
    # Save results
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        results_file = output_path / "clustering_results.json"
        save_json_safe(results, str(results_file))
        
        if visualize:
            _generate_clustering_visualization(results, routes, output_path)
    
    print(f"✅ Clustering complete! Found {results['optimal_k']} clusters")
    return results


def quick_retrieve(
    query: str,
    metadata_db: str,
    route_db_dir: str,
    embedding_cache: str,
    top_k: int = 10,
    output_dir: Optional[str] = None,
    config_path: Optional[str] = None,
    visualize: bool = True
) -> List[Dict[str, Any]]:
    """
    Quick retrieval of similar strategies.
    
    Args:
        query: Text query for retrieval
        metadata_db: Path to metadata database
        route_db_dir: Directory containing route database
        embedding_cache: Path to embedding cache
        top_k: Number of top results to return
        output_dir: Output directory for results
        config_path: Path to configuration file
        visualize: Whether to generate visualizations
        
    Returns:
        List of retrieval results
    """
    # Load configuration
    config = load_config(config_path) if config_path else Config()
    
    # Initialize retriever
    from .retrieval.retriever import StrategyRetriever
    from .models.sentence_transformer_embedder import SentenceTransformerEmbedder
    
    embedder = SentenceTransformerEmbedder(model_name="all-MiniLM-L6-v2")
    
    retriever = StrategyRetriever(
        metadata_db_path=metadata_db,
        route_db_dir=route_db_dir,
        embedding_cache_path=embedding_cache,
        embedder=embedder
    )
    
    print(f"🔍 Retrieving strategies for: '{query}'")
    
    # Create query structure
    query_data = {
        "queries": [{
            "query": {
                "natural_language_description": query
            }
        }]
    }
    
    # Perform retrieval
    results = retriever.retrieve_complex(query_data, top_k=top_k)
    
    # Save results
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        results_file = output_path / "retrieval_results.json"
        save_json_safe(results, str(results_file))
        
        if visualize:
            _generate_retrieval_visualization(results, output_path)
    
    print(f"✅ Retrieval complete! Found {len(results)} results")
    return results


def quick_analyze(
    routes_path: Union[str, List[Dict[str, Any]]],
    functions_dir: str,
    output_dir: Optional[str] = None,
    config_path: Optional[str] = None
) -> Dict[str, Any]:
    """
    Quick analysis of synthesis routes (annotation + clustering).
    
    Args:
        routes_path: Path to routes directory/file or list of route dictionaries
        functions_dir: Directory containing strategy functions
        output_dir: Output directory for results
        config_path: Path to configuration file
        
    Returns:
        Analysis results dictionary
    """
    # Load configuration
    config = load_config(config_path) if config_path else Config()
    
    # Initialize API
    api = SynthStrategyAPI()
    
    # Load routes
    if isinstance(routes_path, str):
        if Path(routes_path).is_dir():
            routes = _load_routes_from_directory(routes_path)
        else:
            routes = load_json_safe(routes_path, [])
    else:
        routes = routes_path
    
    if not routes:
        raise ValueError("No routes found to analyze")
    
    print(f"📊 Analyzing {len(routes)} routes...")
    
    # Step 1: Annotate routes
    print("Step 1: Annotating routes...")
    annotated_routes = api.annotate_strategies(
        routes=routes,
        functions_dir=functions_dir,
        num_workers=config.defaults.workers
    )
    
    # Step 2: Cluster strategies
    print("Step 2: Clustering strategies...")
    clustering_results = api.cluster_strategies(
        annotated_routes=annotated_routes,
        code_dir=functions_dir
    )
    
    # Combine results
    analysis_results = {
        "annotation": {
            "total_routes": len(routes),
            "annotated_routes": len(annotated_routes)
        },
        "clustering": clustering_results
    }
    
    # Save results
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        # Save annotated routes
        annotated_file = output_path / "annotated_routes.json"
        save_json_safe(annotated_routes, str(annotated_file))
        
        # Save analysis results
        analysis_file = output_path / "analysis_results.json"
        save_json_safe(analysis_results, str(analysis_file))
        
        # Generate visualizations
        if clustering_results.get('status') == 'success':
            _generate_clustering_visualization(clustering_results, annotated_routes, output_path)
    
    print("✅ Analysis complete!")
    return analysis_results


def quick_compare(
    routes1_path: Union[str, List[Dict[str, Any]]],
    routes2_path: Union[str, List[Dict[str, Any]]],
    functions_dir: str,
    output_dir: Optional[str] = None,
    config_path: Optional[str] = None
) -> Dict[str, Any]:
    """
    Quick comparison of two sets of synthesis routes.
    
    Args:
        routes1_path: Path to first routes directory/file or list of route dictionaries
        routes2_path: Path to second routes directory/file or list of route dictionaries
        functions_dir: Directory containing strategy functions
        output_dir: Output directory for results
        config_path: Path to configuration file
        
    Returns:
        Comparison results dictionary
    """
    print("🔄 Comparing two route sets...")
    
    # Analyze both sets
    results1 = quick_analyze(routes1_path, functions_dir, config_path=config_path)
    results2 = quick_analyze(routes2_path, functions_dir, config_path=config_path)
    
    # Compare results
    comparison = {
        "set1": results1,
        "set2": results2,
        "comparison": {
            "route_count_diff": results1["annotation"]["total_routes"] - results2["annotation"]["total_routes"],
            "solved_route_diff": results1["annotation"]["solved_routes"] - results2["annotation"]["solved_routes"],
            "cluster_count_diff": (
                results1["clustering"].get("optimal_k", 0) - 
                results2["clustering"].get("optimal_k", 0)
            )
        }
    }
    
    # Save results
    if output_dir:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        
        comparison_file = output_path / "comparison_results.json"
        save_json_safe(comparison, str(comparison_file))
    
    print("✅ Comparison complete!")
    return comparison


def _load_routes_from_directory(directory: str) -> List[Dict[str, Any]]:
    """Load all routes from a directory of JSON files."""
    routes = []
    json_files = find_json_files(directory)
    
    for file_path in json_files:
        file_routes = load_json_safe(file_path, [])
        if isinstance(file_routes, list):
            routes.extend(file_routes)
        else:
            routes.append(file_routes)
    
    return routes


def _generate_clustering_visualization(
    results: Dict[str, Any], 
    routes: List[Dict[str, Any]], 
    output_dir: Path
) -> None:
    """Generate clustering visualizations."""
    try:
        import subprocess
        import tempfile
        
        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            import json
            json.dump(results, f)
            results_file = f.name
        
        routes_dir = None
        if routes:
            routes_dir = tempfile.mkdtemp()
            with open(Path(routes_dir) / "routes.json", 'w') as f:
                json.dump(routes, f)
        
        # Run visualization script
        cmd = [
            "python", 
            "scripts/visualise_clustering_results.py",
            "--analysis_file", results_file,
            "--output_dir", str(output_dir)
        ]
        
        if routes_dir:
            cmd.extend(["--annotated_dir", routes_dir])
        
        subprocess.run(cmd, check=True)
        print(f"🎨 Clustering visualizations saved to: {output_dir}")
        
    except Exception as e:
        print(f"⚠️  Could not generate clustering visualizations: {e}")


def _generate_retrieval_visualization(
    results: List[Dict[str, Any]], 
    output_dir: Path
) -> None:
    """Generate retrieval visualizations."""
    try:
        import subprocess
        import tempfile
        
        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            import json
            json.dump(results, f)
            results_file = f.name
        
        # Run visualization script
        cmd = [
            "python",
            "scripts/visualise_retrieval_results.py",
            "--retrieval_file", results_file,
            "--output_dir", str(output_dir)
        ]
        
        subprocess.run(cmd, check=True)
        print(f"🎨 Retrieval visualizations saved to: {output_dir}")
        
    except Exception as e:
        print(f"⚠️  Could not generate retrieval visualizations: {e}")


# Convenience aliases
cluster = quick_cluster
retrieve = quick_retrieve
analyze = quick_analyze
compare = quick_compare
