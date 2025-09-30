"""
Interactive notebook interface for SynthStrategy.

This module provides Jupyter notebook-friendly functions for interactive
exploration of clustering and retrieval features.
"""

import json
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import List, Dict, Any, Optional, Union
import ipywidgets as widgets
from IPython.display import display, HTML, Image
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

from .api import SynthStrategyAPI
from .quick import quick_cluster, quick_retrieve, quick_analyze
from .config import load_config, Config
from .utils import load_json_safe, find_json_files


class InteractiveSynthStrategy:
    """
    Interactive interface for SynthStrategy in Jupyter notebooks.
    """
    
    def __init__(self, config_path: Optional[str] = None):
        """
        Initialize interactive interface.
        
        Args:
            config_path: Path to configuration file
        """
        self.config = load_config(config_path) if config_path else Config()
        self.api = SynthStrategyAPI()
        self.routes = None
        self.clustering_results = None
        self.retrieval_results = None
        
    def load_routes(self, routes_path: Union[str, List[Dict[str, Any]]]) -> None:
        """
        Load routes for analysis.
        
        Args:
            routes_path: Path to routes directory/file or list of route dictionaries
        """
        if isinstance(routes_path, str):
            if Path(routes_path).is_dir():
                self.routes = self._load_routes_from_directory(routes_path)
            else:
                self.routes = load_json_safe(routes_path, [])
        else:
            self.routes = routes_path
        
        print(f"📊 Loaded {len(self.routes)} routes")
        self._display_route_summary()
    
    def interactive_cluster(
        self, 
        functions_dir: str,
        max_clusters: int = 15,
        silhouette_threshold: float = 0.3
    ) -> None:
        """
        Interactive clustering with parameter tuning.
        
        Args:
            functions_dir: Directory containing strategy functions
            max_clusters: Maximum number of clusters to consider
            silhouette_threshold: Minimum silhouette score threshold
        """
        if not self.routes:
            raise ValueError("No routes loaded. Call load_routes() first.")
        
        print("🔬 Starting interactive clustering...")
        
        # Create parameter widgets
        max_clusters_widget = widgets.IntSlider(
            value=max_clusters,
            min=2,
            max=min(20, len(self.routes) // 2),
            step=1,
            description='Max Clusters:',
            style={'description_width': 'initial'}
        )
        
        threshold_widget = widgets.FloatSlider(
            value=silhouette_threshold,
            min=0.0,
            max=1.0,
            step=0.05,
            description='Silhouette Threshold:',
            style={'description_width': 'initial'}
        )
        
        cluster_button = widgets.Button(
            description='Run Clustering',
            button_style='success',
            icon='play'
        )
        
        output = widgets.Output()
        
        def on_cluster_click(b):
            with output:
                output.clear_output()
                print("🔄 Running clustering...")
                
                # Update config
                self.config.clustering.max_clusters = max_clusters_widget.value
                self.config.clustering.silhouette_threshold = threshold_widget.value
                
                # Run clustering
                self.clustering_results = self.api.cluster_strategies(
                    annotated_routes=self.routes,
                    code_dir=functions_dir
                )
                
                if self.clustering_results.get('status') == 'success':
                    print(f"✅ Clustering complete! Found {self.clustering_results['optimal_k']} clusters")
                    self._display_clustering_results()
                else:
                    print(f"❌ Clustering failed: {self.clustering_results.get('reason')}")
        
        cluster_button.on_click(on_cluster_click)
        
        # Display widgets
        display(widgets.VBox([
            widgets.HBox([max_clusters_widget, threshold_widget]),
            cluster_button,
            output
        ]))
    
    def interactive_retrieve(
        self,
        metadata_db: str,
        route_db_dir: str,
        embedding_cache: str,
        top_k: int = 10
    ) -> None:
        """
        Interactive strategy retrieval.
        
        Args:
            metadata_db: Path to metadata database
            route_db_dir: Directory containing route database
            embedding_cache: Path to embedding cache
            top_k: Number of top results to return
        """
        print("🔍 Starting interactive retrieval...")
        
        # Create query widget
        query_widget = widgets.Textarea(
            value="oxidation strategy",
            placeholder="Enter your strategy query...",
            description='Query:',
            style={'description_width': 'initial'},
            layout=widgets.Layout(width='100%', height='100px')
        )
        
        top_k_widget = widgets.IntSlider(
            value=top_k,
            min=1,
            max=50,
            step=1,
            description='Top K:',
            style={'description_width': 'initial'}
        )
        
        retrieve_button = widgets.Button(
            description='Search',
            button_style='info',
            icon='search'
        )
        
        output = widgets.Output()
        
        def on_retrieve_click(b):
            with output:
                output.clear_output()
                print("🔄 Searching...")
                
                # Run retrieval
                self.retrieval_results = quick_retrieve(
                    query=query_widget.value,
                    metadata_db=metadata_db,
                    route_db_dir=route_db_dir,
                    embedding_cache=embedding_cache,
                    top_k=top_k_widget.value,
                    visualize=False
                )
                
                print(f"✅ Found {len(self.retrieval_results)} results")
                self._display_retrieval_results()
        
        retrieve_button.on_click(on_retrieve_click)
        
        # Display widgets
        display(widgets.VBox([
            query_widget,
            top_k_widget,
            retrieve_button,
            output
        ]))
    
    def visualize_clusters(self) -> None:
        """Visualize clustering results."""
        if not self.clustering_results or self.clustering_results.get('status') != 'success':
            print("❌ No clustering results available. Run clustering first.")
            return
        
        # Create cluster size bar chart
        cluster_sizes = {}
        for route_idx, cluster_id in self.clustering_results['cluster_allocations'].items():
            cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Cluster sizes
        clusters = sorted(cluster_sizes.keys())
        sizes = [cluster_sizes[c] for c in clusters]
        
        ax1.bar(clusters, sizes)
        ax1.set_xlabel('Cluster ID')
        ax1.set_ylabel('Number of Routes')
        ax1.set_title('Cluster Sizes')
        ax1.grid(True, alpha=0.3)
        
        # Top defining features
        top_features = []
        for cluster_id in clusters:
            features = self.clustering_results['cluster_defining_features'].get(str(cluster_id), [])
            if features:
                top_features.append(features[0][0])  # First feature name
        
        ax2.bar(range(len(top_features)), [1] * len(top_features))
        ax2.set_xlabel('Cluster')
        ax2.set_ylabel('Top Feature')
        ax2.set_title('Top Defining Features by Cluster')
        ax2.set_xticks(range(len(top_features)))
        ax2.set_xticklabels(top_features, rotation=45, ha='right')
        
        plt.tight_layout()
        plt.show()
    
    def visualize_retrieval(self) -> None:
        """Visualize retrieval results."""
        if not self.retrieval_results:
            print("❌ No retrieval results available. Run retrieval first.")
            return
        
        # Create similarity score visualization
        scores = [result.get('similarity_score', 0) for result in self.retrieval_results]
        route_ids = [result.get('route_id', f'Route {i}') for i, result in enumerate(self.retrieval_results)]
        
        fig, ax = plt.subplots(figsize=(12, 6))
        bars = ax.bar(range(len(scores)), scores)
        ax.set_xlabel('Retrieved Routes')
        ax.set_ylabel('Similarity Score')
        ax.set_title('Retrieval Results - Similarity Scores')
        ax.set_xticks(range(len(route_ids)))
        ax.set_xticklabels(route_ids, rotation=45, ha='right')
        
        # Color bars by score
        for i, (bar, score) in enumerate(zip(bars, scores)):
            if score > 0.8:
                bar.set_color('green')
            elif score > 0.6:
                bar.set_color('orange')
            else:
                bar.set_color('red')
        
        plt.tight_layout()
        plt.show()
    
    def explore_cluster(self, cluster_id: int) -> None:
        """
        Explore a specific cluster in detail.
        
        Args:
            cluster_id: ID of cluster to explore
        """
        if not self.clustering_results or self.clustering_results.get('status') != 'success':
            print("❌ No clustering results available. Run clustering first.")
            return
        
        cluster_id_str = str(cluster_id)
        if cluster_id_str not in self.clustering_results['cluster_defining_features']:
            print(f"❌ Cluster {cluster_id} not found.")
            return
        
        print(f"🔍 Exploring Cluster {cluster_id}")
        print("=" * 50)
        
        # Get cluster routes
        cluster_routes = []
        for route_idx, c_id in self.clustering_results['cluster_allocations'].items():
            if c_id == cluster_id:
                cluster_routes.append(int(route_idx))
        
        print(f"📊 Cluster {cluster_id} contains {len(cluster_routes)} routes")
        
        # Display defining features
        features = self.clustering_results['cluster_defining_features'][cluster_id_str]
        print(f"\n🎯 Top 5 Defining Features:")
        for i, (feature_name, score) in enumerate(features[:5], 1):
            print(f"  {i}. {feature_name} (score: {score:.4f})")
            
            # Show docstring if available
            docstring = self.clustering_results['function_docstrings'].get(feature_name, "")
            if docstring:
                print(f"     {docstring[:100]}...")
        
        # Show sample routes
        print(f"\n📋 Sample Routes in Cluster {cluster_id}:")
        for i, route_idx in enumerate(cluster_routes[:5]):
            if route_idx < len(self.routes):
                route = self.routes[route_idx]
                print(f"  {i+1}. Route {route_idx}: {route.get('metadata', {}).get('target_smiles', 'N/A')[:50]}...")
    
    def _display_route_summary(self) -> None:
        """Display summary of loaded routes."""
        if not self.routes:
            return
        
        summary_html = f"""
        <div style="background-color: #f0f0f0; padding: 10px; border-radius: 5px; margin: 10px 0;">
            <h3>📊 Route Summary</h3>
            <p><strong>Total Routes:</strong> {len(self.routes)}</p>
        </div>
        """
        display(HTML(summary_html))
    
    def _display_clustering_results(self) -> None:
        """Display clustering results."""
        if not self.clustering_results:
            return
        
        optimal_k = self.clustering_results['optimal_k']
        cluster_allocations = self.clustering_results['cluster_allocations']
        
        # Count routes per cluster
        cluster_sizes = {}
        for route_idx, cluster_id in cluster_allocations.items():
            cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
        
        results_html = f"""
        <div style="background-color: #e8f5e8; padding: 10px; border-radius: 5px; margin: 10px 0;">
            <h3>✅ Clustering Results</h3>
            <p><strong>Optimal Clusters:</strong> {optimal_k}</p>
            <p><strong>Total Routes:</strong> {len(cluster_allocations)}</p>
            <p><strong>Average Cluster Size:</strong> {len(cluster_allocations)/optimal_k:.1f}</p>
        </div>
        """
        display(HTML(results_html))
        
        # Show cluster sizes
        print("\n📊 Cluster Sizes:")
        for cluster_id in sorted(cluster_sizes.keys()):
            print(f"  Cluster {cluster_id}: {cluster_sizes[cluster_id]} routes")
    
    def _display_retrieval_results(self) -> None:
        """Display retrieval results."""
        if not self.retrieval_results:
            return
        
        results_html = f"""
        <div style="background-color: #e8f0ff; padding: 10px; border-radius: 5px; margin: 10px 0;">
            <h3>🔍 Retrieval Results</h3>
            <p><strong>Found:</strong> {len(self.retrieval_results)} results</p>
        </div>
        """
        display(HTML(results_html))
        
        # Show top results
        print("\n📋 Top Results:")
        for i, result in enumerate(self.retrieval_results[:5], 1):
            route_id = result.get('route_id', f'Result {i}')
            score = result.get('similarity_score', 0)
            print(f"  {i}. {route_id} (score: {score:.3f})")
    
    def _load_routes_from_directory(self, directory: str) -> List[Dict[str, Any]]:
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


# Convenience function for quick setup
def setup_interactive(config_path: Optional[str] = None) -> InteractiveSynthStrategy:
    """
    Quick setup for interactive SynthStrategy in notebooks.
    
    Args:
        config_path: Path to configuration file
        
    Returns:
        InteractiveSynthStrategy instance
    """
    return InteractiveSynthStrategy(config_path)
