# src/steerable_retro/evaluate.py

import argparse
import asyncio
import json
import pickle

from pathlib import Path
from typing import Dict, List, Any
from tqdm import tqdm

from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np

# Local imports
from .retrieval.llm_analysis import generate_strategy_description, rewrite_query
from .retrieval.retriever import StrategyRetriever
from .models import Embedder, OpenAIEmbedder, SciBertEmbedder, SentenceTransformerEmbedder

# --- Helper Functions ---

def load_json(path: Path) -> Any:
    """Loads a JSON file."""
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except (FileNotFoundError, json.JSONDecodeError):
        return {} if path.suffix == '.json' else []

def save_json(data: Any, path: Path):
    """Saves data to a JSON file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(data, f, indent=2)

def get_embedder(model_name: str) -> Embedder:
    """Factory function to get an embedder instance by model name."""
    if 'text-embedding' in model_name:
        from .models.openai_embedder import OpenAIEmbedder
        return OpenAIEmbedder(model_name=model_name)
    elif 'scibert' in model_name:
        from .models.scibert_embedder import SciBertEmbedder
        return SciBertEmbedder(model_name=model_name)
    else: # Assume sentence-transformer
        from .models.sentence_transformer_embedder import SentenceTransformerEmbedder
        return SentenceTransformerEmbedder(model_name=model_name)

def precompute_corpus_embeddings(embedder: Embedder, metadata_path: Path, cache_path: Path):
    """Generates and caches embeddings for the function metadata corpus if not already present."""
    if cache_path.exists():
        print(f"Corpus embeddings already exist at '{cache_path}'. Skipping generation.")
        return

    print(f"Corpus embeddings not found. Generating for model '{embedder.model_name}'...")
    metadata = load_json(metadata_path)
    if not metadata:
        raise ValueError(f"Could not load metadata from {metadata_path}")

    descriptions = [item.get('description', '') for item in metadata]
    
    print(f"Generating embeddings for {len(descriptions)} function descriptions...")
    embeddings = embedder.get_embeddings(descriptions) # Batch processing

    print(f"Saving embeddings to '{cache_path}'...")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'wb') as f:
        pickle.dump(embeddings, f)
    print("Corpus embeddings generated and saved successfully.")

def run_description_generation_stage(benchmark: List[Dict], model: str, cache_path: Path, force: bool, max_workers: int = None) -> Dict[str, str]:
    """Stage 1: Generate natural language descriptions for each route in the benchmark."""
    print("\n--- Stage 1: Generating Strategy Descriptions ---")
    cached_descriptions = {} if force else load_json(cache_path)
    
    items_to_process = []
    for item in benchmark:
        route_id = item['selected_route_id']
        # if route_id not in cached_descriptions:
        #     items_to_process.append(item)
    
    if not items_to_process:
        print("All descriptions found in cache.")
        return cached_descriptions
    
    print(f"Generating descriptions for {len(items_to_process)} new routes...")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Submit all tasks
        future_to_item = {
            executor.submit(generate_strategy_description, item['route_data'], model): item 
            for item in items_to_process
        }
        
        # Process completed tasks with progress bar
        for future in tqdm(as_completed(future_to_item), total=len(future_to_item), desc="Generating descriptions"):
            item = future_to_item[future]
            try:
                description = future.result()
                cached_descriptions[item['selected_route_id']] = description
            except Exception as e:
                print(f"Error generating description for route {item['selected_route_id']}: {e}")
    
    save_json(cached_descriptions, cache_path)
    print(f"Descriptions saved to '{cache_path}'.")
    return cached_descriptions

def run_query_rewriting_stage(descriptions: Dict[str, str], model: str, cache_path: Path, force: bool, max_workers: int = None) -> Dict[str, Dict]:
    """Stage 2: Rewrite natural language descriptions into structured JSON queries."""
    print("\n--- Stage 2: Rewriting Descriptions to Queries ---")
    cached_queries = {} if force else load_json(cache_path)
    
    items_to_process = []
    for route_id, desc in descriptions.items():
        if route_id not in cached_queries and desc:
            items_to_process.append((route_id, desc))
    
    if not items_to_process:
        print("All queries found in cache.")
        return cached_queries
    
    print(f"Rewriting {len(items_to_process)} new descriptions into queries...")
    
    with ThreadPoolExecutor(max_workers=1) as executor:
        # Submit all tasks
        future_to_route = {
            executor.submit(rewrite_query, desc, model): route_id 
            for route_id, desc in items_to_process
        }
        
        # Process completed tasks with progress bar
        for future in tqdm(as_completed(future_to_route), total=len(future_to_route), desc="Rewriting queries"):
            route_id = future_to_route[future]
            try:
                query = future.result()
                if query:  # Only save if the LLM returned a valid query
                    cached_queries[route_id] = query
            except Exception as e:
                print(f"Error rewriting query for route {route_id}: {e}")
    
    save_json(cached_queries, cache_path)
    print(f"Rewritten queries saved to '{cache_path}'.")
    return cached_queries

def run_retrieval_stage(retriever: StrategyRetriever, queries: Dict[str, Dict], top_k: int) -> Dict[str, List[str]]:
    """Stage 3: Run retrieval for each query and get top-k results."""
    print("\n--- Stage 3: Running Retrieval ---")
    retrieval_results = {}
    for route_id, query in queries.items():
        print(f"Retrieving for route_id: {route_id}")
        if not query or 'queries' not in query:
            print(f"  -> Skipping due to invalid query.")
            retrieval_results[route_id] = []
            continue
        
        results = retriever.retrieve_complex(query, top_k=top_k, top_n_functions=25)
        retrieved_ids = [res['route_id'] for res in results]
        retrieval_results[route_id] = retrieved_ids
        print(f"  -> Found {len(retrieved_ids)} results.")
    return retrieval_results

def calculate_and_show_accuracy(results: Dict[str, List[str]], top_k_values: List[int]):
    """Stage 4: Calculate and display Top-K accuracy."""
    print("\n--- Stage 4: Calculating Accuracy ---")
    total_queries = len(results)
    if total_queries == 0:
        print("No results to evaluate.")
        return

    hits = {k: 0 for k in top_k_values}
    
    for ground_truth_id, retrieved_ids in results.items():
        for k in top_k_values:
            if ground_truth_id in retrieved_ids[:k]:
                hits[k] += 1
    
    print("\n--- Evaluation Summary ---")
    print(f"Total Routes Evaluated: {total_queries}")
    for k in sorted(hits.keys()):
        accuracy = (hits[k] / total_queries) * 100 if total_queries > 0 else 0
        print(f"  Top-{k} Accuracy: {hits[k]}/{total_queries} ({accuracy:.2f}%)")
    print("--------------------------")


async def main():
    """Main function to orchestrate the evaluation pipeline."""
    parser = argparse.ArgumentParser(description="End-to-end strategy retrieval evaluation pipeline.")
    
    # --- Path Arguments ---
    parser.add_argument("--benchmark_path", type=str, default="data/benchmark/merged_bench.json", help="Path to the benchmark JSON file.")
    parser.add_argument("--metadata_path", type=str, default="data/retrieval/function_metadata_database.json", help="Path to function metadata.")
    parser.add_argument("--route_db_dir", type=str, default="data/benchmark/paroutesn5_structured", help="Directory containing route data.")
    parser.add_argument("--output_dir", type=str, default="data/benchmark/eval_runs", help="Directory to save caches and final results.")

    # --- Model Arguments ---
    parser.add_argument("--embedding_model", type=str, default="text-embedding-3-large", help="Name of the embedding model to use.")
    parser.add_argument("--description_model", type=str, default="google/gemini-2.5-pro", help="LLM for generating strategy descriptions.")
    parser.add_argument("--rewriter_model", type=str, default="google/gemini-2.5-pro", help="LLM for rewriting descriptions to queries.")

    # --- Evaluation Arguments ---
    parser.add_argument("--top_k", nargs='+', type=int, default=[1, 2, 3, 4, 5, 10, 15, 20, 30, 40, 50], help="List of K values for Top-K accuracy.")
    
    # --- Control Arguments ---
    parser.add_argument("--force_description", action="store_true", help="Force regeneration of all strategy descriptions.")
    parser.add_argument("--force_rewrite", action="store_true", help="Force regeneration of all rewritten queries.")

    args = parser.parse_args()

    # --- Setup Paths ---
    base_path = Path(__file__).resolve().parents[2] # Root of steerable_retro project
    benchmark_path = base_path / args.benchmark_path
    metadata_path = base_path / args.metadata_path
    route_db_dir = base_path / args.route_db_dir
    output_dir = base_path / args.output_dir
    output_dir.mkdir(exist_ok=True)

    corpus_embedding_cache_dir = base_path / "data/retrieval_embeddings"
    corpus_embedding_cache_path = corpus_embedding_cache_dir / f"{args.embedding_model.replace('/', '_')}.pck"
    
    description_cache_path = output_dir / f"merged_descriptions.json"
    rewritten_query_cache_path = output_dir / f"queries/queries_{args.rewriter_model.replace('/', '_')}_2.json"
    final_results_path = output_dir / f"results/retrieval_results_{args.embedding_model.replace('/', '_')}.json"

    # --- Run Pipeline ---
    benchmark_data = load_json(benchmark_path)
    if not benchmark_data:
        print(f"Error: Benchmark file not found or empty at '{benchmark_path}'")
        return

    # 0. Pre-computation and Initialization
    embedder = get_embedder(args.embedding_model)
    precompute_corpus_embeddings(embedder, metadata_path, corpus_embedding_cache_path)
    
    retriever = StrategyRetriever(
        metadata_db_path=str(metadata_path),
        route_db_dir=str(route_db_dir),
        embedding_cache_path=str(corpus_embedding_cache_path),
        embedder=embedder
    )

    # 1. Generate Descriptions
    descriptions = run_description_generation_stage(benchmark_data, args.description_model, description_cache_path, args.force_description)
    # print(f"Generated {descriptions} strategy descriptions.")
    # 2. Rewrite Queries
    queries = run_query_rewriting_stage(descriptions, args.rewriter_model, rewritten_query_cache_path, args.force_rewrite)

    # 3. Run Retrieval
    max_k = max(args.top_k)

    
    retrieval_results = run_retrieval_stage(retriever, queries, top_k=max_k)
    save_json(retrieval_results, final_results_path)
    print(f"\nFinal retrieval results saved to '{final_results_path}'")

    # 4. Calculate Accuracy
    calculate_and_show_accuracy(retrieval_results, args.top_k)


if __name__ == "__main__":
    # Make sure you have set your OPENAI_API_KEY (or other necessary keys for litellm)
    # as an environment variable before running.
    asyncio.run(main())