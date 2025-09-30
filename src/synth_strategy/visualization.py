"""
Integrated visualization for SynthStrategy.

This module provides built-in visualization capabilities for clustering and retrieval results,
integrated directly into the main workflow.
"""

import json
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd

from .utils import save_json_safe, ensure_directory


class ClusteringVisualizer:
    """Visualizer for clustering results."""
    
    def __init__(self, results: Dict[str, Any], routes: List[Dict[str, Any]] = None):
        """
        Initialize clustering visualizer.
        
        Args:
            results: Clustering results dictionary
            routes: List of routes used for clustering
        """
        self.results = results
        self.routes = routes or []
        
    def visualize(self, output_dir: Optional[str] = None, interactive: bool = False) -> None:
        """
        Generate comprehensive clustering visualizations.
        
        Args:
            output_dir: Directory to save visualizations
            interactive: Whether to create interactive plots
        """
        if self.results.get('status') != 'success':
            print(f"❌ Cannot visualize failed clustering: {self.results.get('reason')}")
            return
        
        if interactive:
            self._create_interactive_plots(output_dir)
        else:
            self._create_static_plots(output_dir)
    
    def _create_static_plots(self, output_dir: Optional[str] = None) -> None:
        """Create static matplotlib plots."""
        # Set style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # Create figure with subplots
        fig = plt.figure(figsize=(20, 12))
        
        # 1. Cluster size distribution
        ax1 = plt.subplot(2, 3, 1)
        self._plot_cluster_sizes(ax1)
        
        # 2. Silhouette scores
        ax2 = plt.subplot(2, 3, 2)
        self._plot_silhouette_scores(ax2)
        
        # 3. Feature importance heatmap
        ax3 = plt.subplot(2, 3, 3)
        self._plot_feature_heatmap(ax3)
        
        # 4. Cluster composition
        ax4 = plt.subplot(2, 3, 4)
        self._plot_cluster_composition(ax4)
        
        # 5. Route distribution
        ax5 = plt.subplot(2, 3, 5)
        self._plot_route_distribution(ax5)
        
        # 6. Summary statistics
        ax6 = plt.subplot(2, 3, 6)
        self._plot_summary_stats(ax6)
        
        plt.tight_layout()
        
        if output_dir:
            output_path = ensure_directory(output_dir)
            plt.savefig(output_path / "clustering_analysis.png", dpi=300, bbox_inches='tight')
            print(f"📊 Clustering visualizations saved to: {output_path}")
        
        plt.show()
    
    def _create_interactive_plots(self, output_dir: Optional[str] = None) -> None:
        """Create interactive plotly plots."""
        # Create subplots
        fig = make_subplots(
            rows=2, cols=3,
            subplot_titles=[
                "Cluster Sizes", "Silhouette Scores", "Feature Importance",
                "Cluster Composition", "Route Distribution", "Summary Statistics"
            ],
            specs=[[{"type": "bar"}, {"type": "bar"}, {"type": "heatmap"}],
                   [{"type": "pie"}, {"type": "scatter"}, {"type": "table"}]]
        )
        
        # Add plots
        self._add_cluster_sizes_plotly(fig, 1, 1)
        self._add_silhouette_scores_plotly(fig, 1, 2)
        self._add_feature_heatmap_plotly(fig, 1, 3)
        self._add_cluster_composition_plotly(fig, 2, 1)
        self._add_route_distribution_plotly(fig, 2, 2)
        self._add_summary_stats_plotly(fig, 2, 3)
        
        fig.update_layout(
            title="Clustering Analysis Dashboard",
            height=800,
            showlegend=False
        )
        
        if output_dir:
            output_path = ensure_directory(output_dir)
            fig.write_html(str(output_path / "clustering_dashboard.html"))
            print(f"📊 Interactive clustering dashboard saved to: {output_path}")
        
        fig.show()
    
    def _plot_cluster_sizes(self, ax) -> None:
        """Plot cluster size distribution."""
        cluster_sizes = {}
        for route_idx, cluster_id in self.results['cluster_allocations'].items():
            cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
        
        clusters = sorted(cluster_sizes.keys())
        sizes = [cluster_sizes[c] for c in clusters]
        
        bars = ax.bar(clusters, sizes, color=plt.cm.Set3(np.linspace(0, 1, len(clusters))))
        ax.set_xlabel('Cluster ID')
        ax.set_ylabel('Number of Routes')
        ax.set_title('Cluster Size Distribution')
        ax.grid(True, alpha=0.3)
        
        # Add value labels on bars
        for bar, size in zip(bars, sizes):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.1,
                   str(size), ha='center', va='bottom')
    
    def _plot_silhouette_scores(self, ax) -> None:
        """Plot silhouette scores for each cluster."""
        # Calculate silhouette scores (simplified)
        cluster_scores = {}
        for cluster_id in self.results['cluster_allocations'].values():
            if cluster_id not in cluster_scores:
                # Use defining features as proxy for silhouette score
                features = self.results['cluster_defining_features'].get(str(cluster_id), [])
                if features:
                    cluster_scores[cluster_id] = features[0][1]  # First feature score
                else:
                    cluster_scores[cluster_id] = 0.5
        
        clusters = sorted(cluster_scores.keys())
        scores = [cluster_scores[c] for c in clusters]
        
        bars = ax.bar(clusters, scores, color=plt.cm.viridis(np.array(scores)))
        ax.set_xlabel('Cluster ID')
        ax.set_ylabel('Silhouette Score')
        ax.set_title('Cluster Quality Scores')
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 1)
    
    def _plot_feature_heatmap(self, ax) -> None:
        """Plot feature importance heatmap."""
        # Extract top features for each cluster
        cluster_features = {}
        for cluster_id, features in self.results['cluster_defining_features'].items():
            cluster_features[int(cluster_id)] = [f[0] for f in features[:5]]
        
        # Create feature-cluster matrix
        all_features = set()
        for features in cluster_features.values():
            all_features.update(features)
        
        all_features = sorted(list(all_features))
        clusters = sorted(cluster_features.keys())
        
        matrix = np.zeros((len(all_features), len(clusters)))
        for i, feature in enumerate(all_features):
            for j, cluster_id in enumerate(clusters):
                if feature in cluster_features[cluster_id]:
                    # Use feature score as value
                    cluster_features_list = self.results['cluster_defining_features'].get(str(cluster_id), [])
                    for f_name, score in cluster_features_list:
                        if f_name == feature:
                            matrix[i, j] = score
                            break
        
        im = ax.imshow(matrix, cmap='YlOrRd', aspect='auto')
        ax.set_xticks(range(len(clusters)))
        ax.set_xticklabels(clusters)
        ax.set_yticks(range(len(all_features)))
        ax.set_yticklabels(all_features, fontsize=8)
        ax.set_xlabel('Cluster ID')
        ax.set_ylabel('Features')
        ax.set_title('Feature Importance by Cluster')
        
        # Add colorbar
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    
    def _plot_cluster_composition(self, ax) -> None:
        """Plot cluster composition pie chart."""
        cluster_sizes = {}
        for route_idx, cluster_id in self.results['cluster_allocations'].items():
            cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
        
        clusters = sorted(cluster_sizes.keys())
        sizes = [cluster_sizes[c] for c in clusters]
        labels = [f'Cluster {c}' for c in clusters]
        
        wedges, texts, autotexts = ax.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
        ax.set_title('Cluster Composition')
        
        # Make percentage text bold
        for autotext in autotexts:
            autotext.set_fontweight('bold')
    
    def _plot_route_distribution(self, ax) -> None:
        """Plot route distribution scatter plot."""
        if not self.routes:
            ax.text(0.5, 0.5, 'No route data available', ha='center', va='center', transform=ax.transAxes)
            ax.set_title('Route Distribution')
            return
        
        # Extract route properties for visualization
        route_depths = []
        route_scores = []
        cluster_ids = []
        
        for route_idx, cluster_id in self.results['cluster_allocations'].items():
            route_idx = int(route_idx)
            if route_idx < len(self.routes):
                route = self.routes[route_idx]
                # Use route depth as x-axis (simplified)
                depth = len(route.get('reactions', []))
                route_depths.append(depth)
                route_scores.append(route.get('score', 0))
                cluster_ids.append(cluster_id)
        
        if route_depths:
            scatter = ax.scatter(route_depths, route_scores, c=cluster_ids, cmap='tab10', alpha=0.6)
            ax.set_xlabel('Route Depth')
            ax.set_ylabel('Route Score')
            ax.set_title('Route Distribution by Cluster')
            ax.grid(True, alpha=0.3)
            
            # Add colorbar
            plt.colorbar(scatter, ax=ax, label='Cluster ID')
    
    def _plot_summary_stats(self, ax) -> None:
        """Plot summary statistics table."""
        ax.axis('off')
        
        # Calculate statistics
        total_routes = len(self.results['cluster_allocations'])
        optimal_k = self.results['optimal_k']
        avg_cluster_size = total_routes / optimal_k
        
        cluster_sizes = {}
        for route_idx, cluster_id in self.results['cluster_allocations'].items():
            cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
        
        min_size = min(cluster_sizes.values())
        max_size = max(cluster_sizes.values())
        
        # Create table data
        table_data = [
            ['Metric', 'Value'],
            ['Total Routes', str(total_routes)],
            ['Optimal Clusters', str(optimal_k)],
            ['Avg Cluster Size', f'{avg_cluster_size:.1f}'],
            ['Min Cluster Size', str(min_size)],
            ['Max Cluster Size', str(max_size)],
            ['Total Features', str(len(self.results['function_docstrings']))]
        ]
        
        table = ax.table(cellText=table_data[1:], colLabels=table_data[0],
                        cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        
        # Style table
        for i in range(len(table_data)):
            for j in range(2):
                cell = table[(i, j)]
                if i == 0:  # Header
                    cell.set_facecolor('#4CAF50')
                    cell.set_text_props(weight='bold', color='white')
                else:
                    cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
        
        ax.set_title('Clustering Summary Statistics', pad=20)
    
    def _add_cluster_sizes_plotly(self, fig, row, col) -> None:
        """Add cluster sizes plot to plotly figure."""
        cluster_sizes = {}
        for route_idx, cluster_id in self.results['cluster_allocations'].items():
            cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
        
        clusters = sorted(cluster_sizes.keys())
        sizes = [cluster_sizes[c] for c in clusters]
        
        fig.add_trace(
            go.Bar(x=clusters, y=sizes, name="Cluster Sizes"),
            row=row, col=col
        )
    
    def _add_silhouette_scores_plotly(self, fig, row, col) -> None:
        """Add silhouette scores plot to plotly figure."""
        cluster_scores = {}
        for cluster_id in self.results['cluster_allocations'].values():
            if cluster_id not in cluster_scores:
                features = self.results['cluster_defining_features'].get(str(cluster_id), [])
                if features:
                    cluster_scores[cluster_id] = features[0][1]
                else:
                    cluster_scores[cluster_id] = 0.5
        
        clusters = sorted(cluster_scores.keys())
        scores = [cluster_scores[c] for c in clusters]
        
        fig.add_trace(
            go.Bar(x=clusters, y=scores, name="Silhouette Scores"),
            row=row, col=col
        )
    
    def _add_feature_heatmap_plotly(self, fig, row, col) -> None:
        """Add feature heatmap to plotly figure."""
        # Extract top features for each cluster
        cluster_features = {}
        for cluster_id, features in self.results['cluster_defining_features'].items():
            cluster_features[int(cluster_id)] = [f[0] for f in features[:5]]
        
        # Create feature-cluster matrix
        all_features = set()
        for features in cluster_features.values():
            all_features.update(features)
        
        all_features = sorted(list(all_features))
        clusters = sorted(cluster_features.keys())
        
        matrix = np.zeros((len(all_features), len(clusters)))
        for i, feature in enumerate(all_features):
            for j, cluster_id in enumerate(clusters):
                if feature in cluster_features[cluster_id]:
                    cluster_features_list = self.results['cluster_defining_features'].get(str(cluster_id), [])
                    for f_name, score in cluster_features_list:
                        if f_name == feature:
                            matrix[i, j] = score
                            break
        
        fig.add_trace(
            go.Heatmap(z=matrix, x=clusters, y=all_features, colorscale='YlOrRd'),
            row=row, col=col
        )
    
    def _add_cluster_composition_plotly(self, fig, row, col) -> None:
        """Add cluster composition pie chart to plotly figure."""
        cluster_sizes = {}
        for route_idx, cluster_id in self.results['cluster_allocations'].items():
            cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
        
        clusters = sorted(cluster_sizes.keys())
        sizes = [cluster_sizes[c] for c in clusters]
        labels = [f'Cluster {c}' for c in clusters]
        
        fig.add_trace(
            go.Pie(labels=labels, values=sizes, name="Cluster Composition"),
            row=row, col=col
        )
    
    def _add_route_distribution_plotly(self, fig, row, col) -> None:
        """Add route distribution scatter plot to plotly figure."""
        if not self.routes:
            fig.add_annotation(
                text="No route data available",
                xref="x", yref="y",
                x=0.5, y=0.5,
                showarrow=False,
                row=row, col=col
            )
            return
        
        # Extract route properties
        route_depths = []
        route_scores = []
        cluster_ids = []
        
        for route_idx, cluster_id in self.results['cluster_allocations'].items():
            route_idx = int(route_idx)
            if route_idx < len(self.routes):
                route = self.routes[route_idx]
                depth = len(route.get('reactions', []))
                route_depths.append(depth)
                route_scores.append(route.get('score', 0))
                cluster_ids.append(cluster_id)
        
        if route_depths:
            fig.add_trace(
                go.Scatter(x=route_depths, y=route_scores, mode='markers',
                          marker=dict(color=cluster_ids, colorscale='tab10'),
                          name="Route Distribution"),
                row=row, col=col
            )
    
    def _add_summary_stats_plotly(self, fig, row, col) -> None:
        """Add summary statistics table to plotly figure."""
        total_routes = len(self.results['cluster_allocations'])
        optimal_k = self.results['optimal_k']
        avg_cluster_size = total_routes / optimal_k
        
        cluster_sizes = {}
        for route_idx, cluster_id in self.results['cluster_allocations'].items():
            cluster_sizes[cluster_id] = cluster_sizes.get(cluster_id, 0) + 1
        
        min_size = min(cluster_sizes.values())
        max_size = max(cluster_sizes.values())
        
        table_data = [
            ['Total Routes', str(total_routes)],
            ['Optimal Clusters', str(optimal_k)],
            ['Avg Cluster Size', f'{avg_cluster_size:.1f}'],
            ['Min Cluster Size', str(min_size)],
            ['Max Cluster Size', str(max_size)],
            ['Total Features', str(len(self.results['function_docstrings']))]
        ]
        
        fig.add_trace(
            go.Table(
                header=dict(values=['Metric', 'Value'], fill_color='lightblue'),
                cells=dict(values=list(zip(*table_data)), fill_color='white')
            ),
            row=row, col=col
        )


class RetrievalVisualizer:
    """Visualizer for retrieval results."""
    
    def __init__(self, results: List[Dict[str, Any]]):
        """
        Initialize retrieval visualizer.
        
        Args:
            results: List of retrieval results
        """
        self.results = results
    
    def visualize(self, output_dir: Optional[str] = None, interactive: bool = False) -> None:
        """
        Generate comprehensive retrieval visualizations.
        
        Args:
            output_dir: Directory to save visualizations
            interactive: Whether to create interactive plots
        """
        if interactive:
            self._create_interactive_plots(output_dir)
        else:
            self._create_static_plots(output_dir)
    
    def _create_static_plots(self, output_dir: Optional[str] = None) -> None:
        """Create static matplotlib plots."""
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
        
        # 1. Similarity scores
        self._plot_similarity_scores(ax1)
        
        # 2. Score distribution
        self._plot_score_distribution(ax2)
        
        # 3. Top results
        self._plot_top_results(ax3)
        
        # 4. Summary statistics
        self._plot_summary_stats(ax4)
        
        plt.tight_layout()
        
        if output_dir:
            output_path = ensure_directory(output_dir)
            plt.savefig(output_path / "retrieval_analysis.png", dpi=300, bbox_inches='tight')
            print(f"📊 Retrieval visualizations saved to: {output_path}")
        
        plt.show()
    
    def _create_interactive_plots(self, output_dir: Optional[str] = None) -> None:
        """Create interactive plotly plots."""
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=[
                "Similarity Scores", "Score Distribution",
                "Top Results", "Summary Statistics"
            ],
            specs=[[{"type": "bar"}, {"type": "histogram"}],
                   [{"type": "bar"}, {"type": "table"}]]
        )
        
        # Add plots
        self._add_similarity_scores_plotly(fig, 1, 1)
        self._add_score_distribution_plotly(fig, 1, 2)
        self._add_top_results_plotly(fig, 2, 1)
        self._add_summary_stats_plotly(fig, 2, 2)
        
        fig.update_layout(
            title="Retrieval Analysis Dashboard",
            height=600,
            showlegend=False
        )
        
        if output_dir:
            output_path = ensure_directory(output_dir)
            fig.write_html(str(output_path / "retrieval_dashboard.html"))
            print(f"📊 Interactive retrieval dashboard saved to: {output_path}")
        
        fig.show()
    
    def _plot_similarity_scores(self, ax) -> None:
        """Plot similarity scores for retrieved results."""
        scores = [result.get('similarity_score', 0) for result in self.results]
        route_ids = [result.get('route_id', f'Result {i}') for i, result in enumerate(self.results)]
        
        bars = ax.bar(range(len(scores)), scores, color=plt.cm.viridis(np.array(scores)))
        ax.set_xlabel('Retrieved Routes')
        ax.set_ylabel('Similarity Score')
        ax.set_title('Retrieval Results - Similarity Scores')
        ax.set_xticks(range(len(route_ids)))
        ax.set_xticklabels(route_ids, rotation=45, ha='right')
        ax.grid(True, alpha=0.3)
        
        # Color bars by score
        for i, (bar, score) in enumerate(zip(bars, scores)):
            if score > 0.8:
                bar.set_color('green')
            elif score > 0.6:
                bar.set_color('orange')
            else:
                bar.set_color('red')
    
    def _plot_score_distribution(self, ax) -> None:
        """Plot distribution of similarity scores."""
        scores = [result.get('similarity_score', 0) for result in self.results]
        
        ax.hist(scores, bins=20, alpha=0.7, color='skyblue', edgecolor='black')
        ax.set_xlabel('Similarity Score')
        ax.set_ylabel('Frequency')
        ax.set_title('Score Distribution')
        ax.grid(True, alpha=0.3)
        
        # Add statistics
        mean_score = np.mean(scores)
        ax.axvline(mean_score, color='red', linestyle='--', label=f'Mean: {mean_score:.3f}')
        ax.legend()
    
    def _plot_top_results(self, ax) -> None:
        """Plot top N results."""
        top_n = min(10, len(self.results))
        top_results = self.results[:top_n]
        
        scores = [result.get('similarity_score', 0) for result in top_results]
        route_ids = [result.get('route_id', f'Result {i}') for i, result in enumerate(top_results)]
        
        bars = ax.barh(range(len(scores)), scores, color=plt.cm.plasma(np.array(scores)))
        ax.set_xlabel('Similarity Score')
        ax.set_ylabel('Rank')
        ax.set_title(f'Top {top_n} Results')
        ax.set_yticks(range(len(route_ids)))
        ax.set_yticklabels([f'#{i+1}' for i in range(len(route_ids))])
        ax.grid(True, alpha=0.3)
        
        # Add value labels
        for i, (bar, score) in enumerate(zip(bars, scores)):
            ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height()/2,
                   f'{score:.3f}', ha='left', va='center')
    
    def _plot_summary_stats(self, ax) -> None:
        """Plot summary statistics table."""
        ax.axis('off')
        
        # Calculate statistics
        scores = [result.get('similarity_score', 0) for result in self.results]
        mean_score = np.mean(scores)
        std_score = np.std(scores)
        max_score = np.max(scores)
        min_score = np.min(scores)
        
        table_data = [
            ['Metric', 'Value'],
            ['Total Results', str(len(self.results))],
            ['Mean Score', f'{mean_score:.3f}'],
            ['Std Score', f'{std_score:.3f}'],
            ['Max Score', f'{max_score:.3f}'],
            ['Min Score', f'{min_score:.3f}'],
            ['High Quality (>0.8)', str(sum(1 for s in scores if s > 0.8))],
            ['Medium Quality (0.6-0.8)', str(sum(1 for s in scores if 0.6 <= s <= 0.8))],
            ['Low Quality (<0.6)', str(sum(1 for s in scores if s < 0.6))]
        ]
        
        table = ax.table(cellText=table_data[1:], colLabels=table_data[0],
                        cellLoc='center', loc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1.2, 1.5)
        
        # Style table
        for i in range(len(table_data)):
            for j in range(2):
                cell = table[(i, j)]
                if i == 0:  # Header
                    cell.set_facecolor('#4CAF50')
                    cell.set_text_props(weight='bold', color='white')
                else:
                    cell.set_facecolor('#f0f0f0' if i % 2 == 0 else 'white')
        
        ax.set_title('Retrieval Summary Statistics', pad=20)
    
    def _add_similarity_scores_plotly(self, fig, row, col) -> None:
        """Add similarity scores plot to plotly figure."""
        scores = [result.get('similarity_score', 0) for result in self.results]
        route_ids = [result.get('route_id', f'Result {i}') for i, result in enumerate(self.results)]
        
        fig.add_trace(
            go.Bar(x=route_ids, y=scores, name="Similarity Scores"),
            row=row, col=col
        )
    
    def _add_score_distribution_plotly(self, fig, row, col) -> None:
        """Add score distribution plot to plotly figure."""
        scores = [result.get('similarity_score', 0) for result in self.results]
        
        fig.add_trace(
            go.Histogram(x=scores, nbinsx=20, name="Score Distribution"),
            row=row, col=col
        )
    
    def _add_top_results_plotly(self, fig, row, col) -> None:
        """Add top results plot to plotly figure."""
        top_n = min(10, len(self.results))
        top_results = self.results[:top_n]
        
        scores = [result.get('similarity_score', 0) for result in top_results]
        route_ids = [result.get('route_id', f'Result {i}') for i, result in enumerate(top_results)]
        
        fig.add_trace(
            go.Bar(x=scores, y=[f'#{i+1}' for i in range(len(route_ids))],
                   orientation='h', name="Top Results"),
            row=row, col=col
        )
    
    def _add_summary_stats_plotly(self, fig, row, col) -> None:
        """Add summary statistics table to plotly figure."""
        scores = [result.get('similarity_score', 0) for result in self.results]
        mean_score = np.mean(scores)
        std_score = np.std(scores)
        max_score = np.max(scores)
        min_score = np.min(scores)
        
        table_data = [
            ['Total Results', str(len(self.results))],
            ['Mean Score', f'{mean_score:.3f}'],
            ['Std Score', f'{std_score:.3f}'],
            ['Max Score', f'{max_score:.3f}'],
            ['Min Score', f'{min_score:.3f}'],
            ['High Quality (>0.8)', str(sum(1 for s in scores if s > 0.8))],
            ['Medium Quality (0.6-0.8)', str(sum(1 for s in scores if 0.6 <= s <= 0.8))],
            ['Low Quality (<0.6)', str(sum(1 for s in scores if s < 0.6))]
        ]
        
        fig.add_trace(
            go.Table(
                header=dict(values=['Metric', 'Value'], fill_color='lightblue'),
                cells=dict(values=list(zip(*table_data)), fill_color='white')
            ),
            row=row, col=col
        )


def visualize_clustering(results: Dict[str, Any], routes: List[Dict[str, Any]] = None, 
                        output_dir: Optional[str] = None, interactive: bool = False) -> None:
    """
    Quick function to visualize clustering results.
    
    Args:
        results: Clustering results dictionary
        routes: List of routes used for clustering
        output_dir: Directory to save visualizations
        interactive: Whether to create interactive plots
    """
    visualizer = ClusteringVisualizer(results, routes)
    visualizer.visualize(output_dir, interactive)


def visualize_retrieval(results: List[Dict[str, Any]], output_dir: Optional[str] = None, 
                       interactive: bool = False) -> None:
    """
    Quick function to visualize retrieval results.
    
    Args:
        results: List of retrieval results
        output_dir: Directory to save visualizations
        interactive: Whether to create interactive plots
    """
    visualizer = RetrievalVisualizer(results)
    visualizer.visualize(output_dir, interactive)
