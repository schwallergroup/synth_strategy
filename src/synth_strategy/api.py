# src/synth_strategy/api.py

import asyncio
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

from .route_annotation.annotation import annotate_routes as perform_annotation
from .clustering.clustering import perform_strategy_clustering

from .llm.function_generation import LM as StrategyGeneratorLM, load_checker

# For running the end-to-end retrieval and evaluation pipeline
# We assume evaluate.py is refactored to expose its core logic
try:
    # This requires refactoring evaluate.py's main() logic into this new function
    from .evaluate import execute_evaluation_pipeline
except ImportError:
    print("Warning: Could not import from evaluate. Retrieval evaluation will not be available.", file=sys.stderr)
    execute_evaluation_pipeline = None
# --- FIX END ---


class SynthStrategyAPI:
    """
    A unified API for strategy annotation, clustering, function generation,
    and retrieval evaluation in the Steerable Retrosynthesis project.

    This class provides a single entry point to the core functionalities,
    allowing for programmatic control over the pipeline stages.

    Example Usage:
        # Initialize the API
        api = SynthStrategyAPI()

        # 1. Annotate routes
        raw_routes = [...] # Load your raw route data
        annotated = api.annotate_strategies(
            routes=raw_routes,
            functions_dir="path/to/strategy_functions"
        )

        # 2. Cluster the annotated strategies
        clusters = api.cluster_strategies(
            annotated_routes=annotated,
            code_dir="path/to/strategy_functions" # For docstrings
        )
        print(f"Found {clusters.get('optimal_k')} strategy clusters.")

        # 3. Generate a new strategy function from a route
        route_to_generate = raw_routes[0]
        model_cfg = {"model_name": "claude-3-5-sonnet-20240620"}
        dep_paths = {
            "fg_path": "path/to/functional_groups.json",
            "reaction_path": "path/to/reactions.json",
            "ring_path": "path/to/rings.json"
        }
        generation_result = asyncio.run(api.generate_strategy_function(
            route=route_to_generate,
            model_config=model_cfg,
            dependency_paths=dep_paths
        ))
        print(f"Generated code passed validation: {generation_result.get('final_passed')}")

        # 4. Run the full evaluation pipeline
        eval_cfg = {"benchmark_path": "...", "output_dir": "...", ...}
        asyncio.run(api.run_retrieval_evaluation(eval_config=eval_cfg))
    """

    def __init__(self):
        """Initializes the SynthStrategyAPI."""
        print("SynthStrategyAPI initialized. Ensure all dependencies are available.")

    def annotate_strategies(
        self, routes: List[Dict[str, Any]], functions_dir: str, num_workers: Optional[int] = None
    ) -> List[Dict[str, Any]]:
        """
        Annotates a list of synthesis routes by applying a set of strategy functions in parallel.

        Args:
            routes: A list of raw synthesis route dictionaries.
            functions_dir: The path to the directory containing the Python strategy function files.
            num_workers: The number of parallel processes to use. Defaults to the number of CPU cores.

        Returns:
            The list of routes, annotated with 'passing_functions' and 'errored_functions'.
        """
        print(f"Starting strategy annotation for {len(routes)} routes...")
        annotated = perform_annotation(routes, functions_dir, num_workers)
        print("Annotation complete.")
        return annotated

    def cluster_strategies(
        self, annotated_routes: List[Dict[str, Any]], code_dir: str
    ) -> Dict[str, Any]:
        """
        Performs strategy-based clustering on a list of pre-annotated synthesis routes.

        Args:
            annotated_routes: A list of routes that have been processed by `annotate_strategies`.
            code_dir: The path to the source code directory of the strategy functions,
                      used to extract docstrings for context.

        Returns:
            A dictionary containing the clustering analysis, including the optimal number
            of clusters, cluster assignments, and defining features for each cluster.
        """
        print(f"Starting clustering for {len(annotated_routes)} annotated routes...")
        results = perform_strategy_clustering(
            routes=annotated_routes,
            code_dir=code_dir,
            pre_annotated=True  # We assume the input is already annotated
        )
        print("Clustering complete.")
        return results

    async def generate_strategy_function(
        self,
        route: Dict[str, Any],
        model_config: Dict[str, Any],
        dependency_paths: Dict[str, str]
    ) -> Dict[str, Any]:
        """
        Generates a Python strategy function for a single synthesis route using an LLM.

        This is an asynchronous operation.

        Args:
            route: A single synthesis route dictionary.
            model_config: Configuration for the LLM (e.g., {'model_name': 'claude-3-5-sonnet-20240620'}).
            dependency_paths: A dictionary with paths to required data files:
                              'fg_path', 'reaction_path', and 'ring_path'.

        Returns:
            A dictionary with the generation results, including the final code, pass/fail status,
            and a history of self-correction attempts.
        """
        if StrategyGeneratorLM is None or load_checker is None:
            raise NotImplementedError("Function generation dependencies could not be imported from 'llm.function_generation'.")

        print(f"Starting function generation for route...")

        # 1. Load chemical dictionaries and checker utility
        try:
            fg, reactions, rings, checker = load_checker(
                dependency_paths['fg_path'],
                dependency_paths['reaction_path'],
                dependency_paths['ring_path']
            )
        except FileNotFoundError as e:
            raise ValueError(f"Could not load dependency files for function generation: {e}")

        # 2. Initialize the Language Model from sequential.py
        lm_generator = StrategyGeneratorLM(
            prompt="synth_strategy.llm.prompts.strategy_extraction",
            model=model_config.get("model_name", "claude-3-5-sonnet-20240620"),
            fg_dict=fg,
            reaction_dict=reactions,
            ring_dict=rings,
            checker=checker,
        )

        # 3. Run the generation and self-correction loop for the single route
        # The original function returns a dictionary {index: result}, so we extract the result.
        result_dict = await lm_generator.run_single_route(d=route, task="", idx=0)
        print("Function generation complete.")
        return result_dict.get(0, {"error": "Generation failed to produce a result."})

    async def run_retrieval_evaluation(self, eval_config: Dict[str, Any]):
        """
        Runs the full end-to-end strategy retrieval and evaluation pipeline.

        This is an asynchronous operation that orchestrates description generation,
        query rewriting, retrieval across different embedding models, and accuracy calculation.
        Results are saved to files as configured.

        Args:
            eval_config: A dictionary containing all configuration parameters needed by the
                         evaluation pipeline, mirroring the arguments in `evaluate.py`.
                         Required keys include: 'benchmark_path', 'metadata_path',
                         'route_db_dir', 'output_dir', 'description_model',
                         'rewriter_model', 'top_k', 'top_n_functions_sweep', etc.
        """
        if execute_evaluation_pipeline is None:
            raise NotImplementedError("Evaluation pipeline could not be imported from 'evaluate.py'.")

        print("Starting the end-to-end retrieval evaluation pipeline...")
        await execute_evaluation_pipeline(eval_config)
        output_dir = eval_config.get('output_dir', 'the configured output directory')
        print(f"Evaluation pipeline finished. Results are saved in '{output_dir}'.")

    def visualize_clustering(self, results: Dict[str, Any], routes: List[Dict[str, Any]] = None, 
                           output_dir: Optional[str] = None, interactive: bool = False) -> None:
        """
        Visualize clustering results with integrated visualization.

        Args:
            results: Clustering results dictionary
            routes: List of routes used for clustering
            output_dir: Directory to save visualizations
            interactive: Whether to create interactive plots
        """
        from .visualization import visualize_clustering
        visualize_clustering(results, routes, output_dir, interactive)

    def visualize_retrieval(self, results: List[Dict[str, Any]], output_dir: Optional[str] = None, 
                          interactive: bool = False) -> None:
        """
        Visualize retrieval results with integrated visualization.

        Args:
            results: List of retrieval results
            output_dir: Directory to save visualizations
            interactive: Whether to create interactive plots
        """
        from .visualization import visualize_retrieval
        visualize_retrieval(results, output_dir, interactive)