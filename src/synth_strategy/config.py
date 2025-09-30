"""
Configuration management for SynthStrategy.

This module provides configuration loading, validation, and default settings
for all SynthStrategy components.
"""

import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from dataclasses import dataclass, field


@dataclass
class DefaultsConfig:
    """Default configuration settings."""
    functions_dir: str = "data/strategy_functions"
    output_dir: str = "results"
    workers: int = 4
    log_level: str = "INFO"


@dataclass
class ClusteringConfig:
    """Clustering-specific configuration."""
    max_clusters: int = 15
    silhouette_threshold: float = 0.3
    min_routes_per_cluster: int = 2
    random_state: int = 42


@dataclass
class RetrievalConfig:
    """Retrieval-specific configuration."""
    embedding_model: str = "all-MiniLM-L6-v2"
    top_k: int = 10
    top_n_functions: Optional[int] = None
    metadata_db: str = "data/function_metadata_database.json"
    route_db_dir: str = "data/routes"
    embedding_cache: str = "data/retrieval_embeddings"


@dataclass
class VisualizationConfig:
    """Visualization-specific configuration."""
    image_format: str = "png"
    examples_per_cluster: int = 10
    dpi: int = 300
    figure_size: tuple = (12, 8)


@dataclass
class Config:
    """Main configuration class."""
    defaults: DefaultsConfig = field(default_factory=DefaultsConfig)
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    retrieval: RetrievalConfig = field(default_factory=RetrievalConfig)
    visualization: VisualizationConfig = field(default_factory=VisualizationConfig)
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert config to dictionary."""
        return {
            "defaults": {
                "functions_dir": self.defaults.functions_dir,
                "output_dir": self.defaults.output_dir,
                "workers": self.defaults.workers,
                "log_level": self.defaults.log_level,
            },
            "clustering": {
                "max_clusters": self.clustering.max_clusters,
                "silhouette_threshold": self.clustering.silhouette_threshold,
                "min_routes_per_cluster": self.clustering.min_routes_per_cluster,
                "random_state": self.clustering.random_state,
            },
            "retrieval": {
                "embedding_model": self.retrieval.embedding_model,
                "top_k": self.retrieval.top_k,
                "top_n_functions": self.retrieval.top_n_functions,
                "metadata_db": self.retrieval.metadata_db,
                "route_db_dir": self.retrieval.route_db_dir,
                "embedding_cache": self.retrieval.embedding_cache,
            },
            "visualization": {
                "image_format": self.visualization.image_format,
                "examples_per_cluster": self.visualization.examples_per_cluster,
                "dpi": self.visualization.dpi,
                "figure_size": self.visualization.figure_size,
            }
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Config":
        """Create config from dictionary."""
        config = cls()
        
        if "defaults" in data:
            defaults = data["defaults"]
            config.defaults = DefaultsConfig(
                functions_dir=defaults.get("functions_dir", config.defaults.functions_dir),
                output_dir=defaults.get("output_dir", config.defaults.output_dir),
                workers=defaults.get("workers", config.defaults.workers),
                log_level=defaults.get("log_level", config.defaults.log_level),
            )
        
        if "clustering" in data:
            clustering = data["clustering"]
            config.clustering = ClusteringConfig(
                max_clusters=clustering.get("max_clusters", config.clustering.max_clusters),
                silhouette_threshold=clustering.get("silhouette_threshold", config.clustering.silhouette_threshold),
                min_routes_per_cluster=clustering.get("min_routes_per_cluster", config.clustering.min_routes_per_cluster),
                random_state=clustering.get("random_state", config.clustering.random_state),
            )
        
        if "retrieval" in data:
            retrieval = data["retrieval"]
            config.retrieval = RetrievalConfig(
                embedding_model=retrieval.get("embedding_model", config.retrieval.embedding_model),
                top_k=retrieval.get("top_k", config.retrieval.top_k),
                top_n_functions=retrieval.get("top_n_functions", config.retrieval.top_n_functions),
                metadata_db=retrieval.get("metadata_db", config.retrieval.metadata_db),
                route_db_dir=retrieval.get("route_db_dir", config.retrieval.route_db_dir),
                embedding_cache=retrieval.get("embedding_cache", config.retrieval.embedding_cache),
            )
        
        if "visualization" in data:
            viz = data["visualization"]
            config.visualization = VisualizationConfig(
                image_format=viz.get("image_format", config.visualization.image_format),
                examples_per_cluster=viz.get("examples_per_cluster", config.visualization.examples_per_cluster),
                dpi=viz.get("dpi", config.visualization.dpi),
                figure_size=tuple(viz.get("figure_size", config.visualization.figure_size)),
            )
        
        return config


def load_config(config_path: Optional[str] = None) -> Config:
    """
    Load configuration from file or return defaults.
    
    Args:
        config_path: Path to YAML configuration file
        
    Returns:
        Config object
    """
    if config_path is None:
        return Config()
    
    config_file = Path(config_path)
    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_path}")
    
    try:
        with open(config_file, 'r') as f:
            data = yaml.safe_load(f)
        
        return Config.from_dict(data)
    
    except yaml.YAMLError as e:
        raise ValueError(f"Invalid YAML in configuration file: {e}")
    except Exception as e:
        raise ValueError(f"Error loading configuration: {e}")


def save_config(config: Config, config_path: str) -> None:
    """
    Save configuration to file.
    
    Args:
        config: Config object to save
        config_path: Path to save configuration file
    """
    config_file = Path(config_path)
    config_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(config_file, 'w') as f:
        yaml.dump(config.to_dict(), f, default_flow_style=False, indent=2)


def create_default_config(config_path: str) -> None:
    """
    Create a default configuration file.
    
    Args:
        config_path: Path to save the default configuration
    """
    config = Config()
    save_config(config, config_path)
    print(f"Default configuration saved to: {config_path}")


def validate_config(config: Config) -> None:
    """
    Validate configuration values.
    
    Args:
        config: Config object to validate
        
    Raises:
        ValueError: If configuration is invalid
    """
    # Validate defaults
    if config.defaults.workers < 1:
        raise ValueError("workers must be at least 1")
    
    if config.defaults.log_level not in ["DEBUG", "INFO", "WARNING", "ERROR"]:
        raise ValueError("log_level must be one of: DEBUG, INFO, WARNING, ERROR")
    
    # Validate clustering
    if config.clustering.max_clusters < 2:
        raise ValueError("max_clusters must be at least 2")
    
    if not 0 <= config.clustering.silhouette_threshold <= 1:
        raise ValueError("silhouette_threshold must be between 0 and 1")
    
    # Validate retrieval
    if config.retrieval.top_k < 1:
        raise ValueError("top_k must be at least 1")
    
    # Validate visualization
    if config.visualization.examples_per_cluster < 1:
        raise ValueError("examples_per_cluster must be at least 1")
    
    if config.visualization.dpi < 72:
        raise ValueError("dpi must be at least 72")
