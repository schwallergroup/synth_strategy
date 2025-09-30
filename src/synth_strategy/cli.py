#!/usr/bin/env python3
"""
Unified CLI for SynthStrategy clustering and retrieval features.

This module provides a single command-line interface for all major functionality
including clustering, retrieval, and visualization.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, Any, List, Optional

from .api import SynthStrategyAPI
from .config import load_config, Config
from .utils import setup_logging, validate_paths


def annotate_command(args) -> None:
    """Handle annotation command."""
    print("📝 Starting route annotation...")
    
    # Load configuration
    config = load_config(args.config) if args.config else Config()
    
    # Initialize API
    api = SynthStrategyAPI()
    
    # Load routes
    if args.input_dir:
        routes = load_routes_from_directory(args.input_dir)
    elif args.input_file:
        routes = load_routes_from_file(args.input_file)
    else:
        print("Error: Must specify either --input-dir or --input-file")
        sys.exit(1)
    
    print(f"📊 Loaded {len(routes)} routes for annotation")
    
    # Perform annotation
    try:
        annotated_routes = api.annotate_strategies(
            routes=routes,
            functions_dir=args.functions_dir or config.defaults.functions_dir,
            num_workers=args.workers or config.defaults.workers
        )
        
        # Save results
        output_path = Path(args.output or config.defaults.output_dir) / "annotated_routes.json"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(annotated_routes, f, indent=2)
        
        print(f"✅ Annotation complete!")
        print(f"📊 Total routes: {len(annotated_routes)}")
        print(f"📁 Results saved to: {output_path}")
        
        # Optionally cluster after annotation
        if args.cluster:
            print("\n🔬 Running clustering on annotated routes...")
            clustering_results = api.cluster_strategies(
                annotated_routes=annotated_routes,
                code_dir=args.functions_dir or config.defaults.functions_dir
            )
            
            if clustering_results.get('status') == 'success':
                cluster_output_path = output_path.parent / "clustering_results.json"
                with open(cluster_output_path, 'w') as f:
                    json.dump(clustering_results, f, indent=2)
                
                print(f"✅ Clustering complete! Found {clustering_results['optimal_k']} clusters")
                print(f"📁 Clustering results saved to: {cluster_output_path}")
                
                # Generate visualization if requested
                if args.visualize:
                    print("🎨 Generating visualizations...")
                    api.visualize_clustering(
                        results=clustering_results,
                        routes=annotated_routes,
                        output_dir=str(output_path.parent / "visualizations"),
                        interactive=args.interactive
                    )
            else:
                print(f"⚠️  Clustering skipped: {clustering_results.get('reason', 'Unknown reason')}")
        
    except Exception as e:
        print(f"❌ Error during annotation: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def cluster_command(args) -> None:
    """Handle clustering command."""
    print("🔬 Starting strategy clustering analysis...")
    
    # Load configuration
    config = load_config(args.config) if args.config else Config()
    
    # Initialize API
    api = SynthStrategyAPI()
    
    # Load routes
    if args.input_dir:
        routes = load_routes_from_directory(args.input_dir)
    elif args.input_file:
        routes = load_routes_from_file(args.input_file)
    else:
        print("Error: Must specify either --input-dir or --input-file")
        sys.exit(1)
    
    print(f"📊 Loaded {len(routes)} routes for clustering")
    
    # Perform clustering
    try:
        results = api.cluster_strategies(
            annotated_routes=routes,
            code_dir=args.functions_dir or config.defaults.functions_dir
        )
        
        if results.get('status') != 'success':
            print(f"❌ Clustering failed: {results.get('reason', 'Unknown error')}")
            sys.exit(1)
        
        # Save results
        output_path = Path(args.output or config.defaults.output_dir) / "clustering_results.json"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"✅ Clustering complete! Found {results['optimal_k']} clusters")
        print(f"📁 Results saved to: {output_path}")
        
        # Generate visualization if requested
        if args.visualize:
            print("🎨 Generating visualizations...")
            visualize_clustering_results(results, routes, output_path.parent)
        
    except Exception as e:
        print(f"❌ Error during clustering: {e}")
        sys.exit(1)


def retrieve_command(args) -> None:
    """Handle retrieval command."""
    print("🔍 Starting strategy retrieval...")
    
    # Load configuration
    config = load_config(args.config) if args.config else Config()
    
    # Initialize retriever
    from .retrieval.retriever import StrategyRetriever
    from .models.sentence_transformer_embedder import SentenceTransformerEmbedder
    
    # Setup embedder
    embedder = SentenceTransformerEmbedder(model_name="all-MiniLM-L6-v2")
    
    # Initialize retriever
    retriever = StrategyRetriever(
        metadata_db_path=args.metadata_db or config.retrieval.metadata_db,
        route_db_dir=args.route_db or config.retrieval.route_db_dir,
        embedding_cache_path=args.embedding_cache or config.retrieval.embedding_cache,
        embedder=embedder
    )
    
    # Perform retrieval
    try:
        # Parse query
        if args.query_file:
            with open(args.query_file, 'r') as f:
                query_data = json.load(f)
        else:
            # Simple text query
            query_data = {
                "queries": [{
                    "query": {
                        "natural_language_description": args.query
                    }
                }]
            }
        
        results = retriever.retrieve_complex(
            query_data, 
            top_k=args.top_k or config.retrieval.top_k
        )
        
        # Save results
        output_path = Path(args.output or config.defaults.output_dir) / "retrieval_results.json"
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_path, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"✅ Retrieval complete! Found {len(results)} results")
        print(f"📁 Results saved to: {output_path}")
        
        # Generate visualization if requested
        if args.visualize:
            print("🎨 Generating visualizations...")
            visualize_retrieval_results(results, output_path.parent)
        
    except Exception as e:
        print(f"❌ Error during retrieval: {e}")
        sys.exit(1)


def visualize_command(args) -> None:
    """Handle visualization command."""
    print("🎨 Starting visualization...")
    
    if args.type == "clustering":
        if not args.input:
            print("Error: --input required for clustering visualization")
            sys.exit(1)
        
        # Load clustering results
        with open(args.input, 'r') as f:
            results = json.load(f)
        
        # Load routes for visualization
        routes = load_routes_from_directory(args.routes_dir) if args.routes_dir else []
        
        visualize_clustering_results(results, routes, args.output)
        
    elif args.type == "retrieval":
        if not args.input:
            print("Error: --input required for retrieval visualization")
            sys.exit(1)
        
        # Load retrieval results
        with open(args.input, 'r') as f:
            results = json.load(f)
        
        visualize_retrieval_results(results, args.output)
    
    else:
        print(f"Error: Unknown visualization type '{args.type}'")
        sys.exit(1)


def load_routes_from_directory(directory: str) -> List[Dict[str, Any]]:
    """Load all routes from a directory of JSON files."""
    import glob
    
    routes = []
    json_files = glob.glob(str(Path(directory) / "*.json"))
    
    for file_path in json_files:
        try:
            with open(file_path, 'r') as f:
                routes.extend(json.load(f))
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load {file_path}: {e}")
    
    return routes


def load_routes_from_file(file_path: str) -> List[Dict[str, Any]]:
    """Load routes from a single JSON file."""
    with open(file_path, 'r') as f:
        return json.load(f)


def visualize_clustering_results(results: Dict[str, Any], routes: List[Dict[str, Any]], output_dir: Path) -> None:
    """Generate clustering visualizations."""
    try:
        # Import visualization script
        import subprocess
        import tempfile
        
        # Create temporary files for the visualization script
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(results, f)
            results_file = f.name
        
        # Create temporary directory for routes if needed
        routes_dir = None
        if routes:
            routes_dir = tempfile.mkdtemp()
            with open(Path(routes_dir) / "routes.json", 'w') as f:
                json.dump(routes, f)
        
        # Run visualization script
        cmd = [
            sys.executable, 
            "scripts/visualise_clustering_results.py",
            "--analysis_file", results_file,
            "--output_dir", str(output_dir)
        ]
        
        if routes_dir:
            cmd.extend(["--annotated_dir", routes_dir])
        
        subprocess.run(cmd, check=True)
        
        print(f"✅ Clustering visualizations saved to: {output_dir}")
        
    except Exception as e:
        print(f"❌ Error generating clustering visualizations: {e}")


def visualize_retrieval_results(results: List[Dict[str, Any]], output_dir: Path) -> None:
    """Generate retrieval visualizations."""
    try:
        # Import visualization script
        import subprocess
        import tempfile
        
        # Create temporary files for the visualization script
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
            json.dump(results, f)
            results_file = f.name
        
        # Run visualization script
        cmd = [
            sys.executable,
            "scripts/visualise_retrieval_results.py",
            "--retrieval_file", results_file,
            "--output_dir", str(output_dir)
        ]
        
        subprocess.run(cmd, check=True)
        
        print(f"✅ Retrieval visualizations saved to: {output_dir}")
        
    except Exception as e:
        print(f"❌ Error generating retrieval visualizations: {e}")


def main():
    """Main CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Unified CLI for SynthStrategy clustering and retrieval",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Annotate routes with strategy functions
  synth-strategy annotate --input-dir data/routes/ --functions-dir data/functions/
  
  # Annotate and cluster in one step
  synth-strategy annotate --input-dir data/routes/ --functions-dir data/functions/ --cluster --visualize
  
  # Cluster already-annotated routes
  synth-strategy cluster --input-dir data/annotated_routes/ --functions-dir data/functions/
  
  # Retrieve strategies with a text query
  synth-strategy retrieve --query "oxidation strategy" --top-k 10
  
  # Visualize clustering results
  synth-strategy visualize --type clustering --input results/clustering_results.json
  
  # Use configuration file
  synth-strategy annotate --config config.yaml --input-dir data/routes/ --cluster
        """
    )
    
    # Global options
    parser.add_argument("--config", "-c", help="Configuration file path")
    parser.add_argument("--verbose", "-v", action="store_true", help="Verbose output")
    
    # Subcommands
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    
    # Annotation command
    annotate_parser = subparsers.add_parser("annotate", help="Annotate routes with strategy functions")
    annotate_parser.add_argument("--input-dir", help="Directory containing route JSON files")
    annotate_parser.add_argument("--input-file", help="Single route JSON file")
    annotate_parser.add_argument("--functions-dir", help="Directory containing strategy functions")
    annotate_parser.add_argument("--output", "-o", help="Output directory")
    annotate_parser.add_argument("--workers", type=int, help="Number of parallel workers")
    annotate_parser.add_argument("--cluster", action="store_true", help="Run clustering after annotation")
    annotate_parser.add_argument("--visualize", action="store_true", help="Generate visualizations (requires --cluster)")
    annotate_parser.add_argument("--interactive", action="store_true", help="Create interactive visualizations")
    
    # Clustering command
    cluster_parser = subparsers.add_parser("cluster", help="Perform strategy clustering")
    cluster_parser.add_argument("--input-dir", help="Directory containing route JSON files")
    cluster_parser.add_argument("--input-file", help="Single route JSON file")
    cluster_parser.add_argument("--functions-dir", help="Directory containing strategy functions")
    cluster_parser.add_argument("--output", "-o", help="Output directory")
    cluster_parser.add_argument("--visualize", action="store_true", help="Generate visualizations")
    
    # Retrieval command
    retrieve_parser = subparsers.add_parser("retrieve", help="Perform strategy retrieval")
    retrieve_parser.add_argument("--query", help="Text query for retrieval")
    retrieve_parser.add_argument("--query-file", help="JSON file containing complex query")
    retrieve_parser.add_argument("--metadata-db", help="Path to metadata database")
    retrieve_parser.add_argument("--route-db", help="Path to route database directory")
    retrieve_parser.add_argument("--embedding-cache", help="Path to embedding cache")
    retrieve_parser.add_argument("--top-k", type=int, help="Number of top results to return")
    retrieve_parser.add_argument("--output", "-o", help="Output directory")
    retrieve_parser.add_argument("--visualize", action="store_true", help="Generate visualizations")
    
    # Visualization command
    viz_parser = subparsers.add_parser("visualize", help="Generate visualizations")
    viz_parser.add_argument("--type", choices=["clustering", "retrieval"], required=True, help="Type of visualization")
    viz_parser.add_argument("--input", "-i", required=True, help="Input results file")
    viz_parser.add_argument("--output", "-o", help="Output directory")
    viz_parser.add_argument("--routes-dir", help="Directory containing routes (for clustering)")
    
    args = parser.parse_args()
    
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    # Setup logging
    setup_logging(verbose=args.verbose)
    
    # Execute command
    if args.command == "annotate":
        annotate_command(args)
    elif args.command == "cluster":
        cluster_command(args)
    elif args.command == "retrieve":
        retrieve_command(args)
    elif args.command == "visualize":
        visualize_command(args)


if __name__ == "__main__":
    main()
