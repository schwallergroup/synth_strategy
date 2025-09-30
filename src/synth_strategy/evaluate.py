import argparse
import asyncio
import csv
import json
import pickle
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Union

import numpy as np
from tqdm import tqdm

from .models import *
from .retrieval.llm_analysis import generate_strategy_description, rewrite_query
from .retrieval.retriever import StrategyRetriever
# NOTE: To make this script runnable standalone, dummy versions of missing classes are provided below.
# In your actual project, you would use your real imports.




# --- CONFIGURATION FOR THE SWEEP ---
HUGGINGFACE_MODELS = [
    "allenai/scibert_scivocab_uncased",
    "bert-base-uncased",
    "FacebookAI/roberta-base",
    "Salesforce/SFR-Embedding-Mistral",
    "FacebookAI/xlm-roberta-base",
    "FacebookAI/xlm-roberta-large",
    "Qwen/Qwen3-Embedding-0.6B",
    "Alibaba-NLP/gte-Qwen2-1.5B-instruct",
    "Qwen/Qwen3-Embedding-4B",
    "intfloat/multilingual-e5-large-instruct",
    "Lajavaness/bilingual-embedding-large",
    "NovaSearch/stella_en_1.5B_v5"
]
EMBEDDING_CONFIGS = [
    # Sentence Transformers
    {"model_name": "all-MiniLM-L6-v2"},
]

def get_config_name(config: dict) -> str:
    """Generates a file-safe name from an embedding configuration dictionary."""
    model_name = config["model_name"].replace('/', '_')
    if "embedding-001" in model_name:
        # For Gemini, include task type and dimension for uniqueness
        task_name = config.get("doc_task_type") or config.get("task_type", "unknown")
        dim = config.get("output_dimensionality", "default")
        return f"{model_name}_{task_name}_{dim}"
    return model_name

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

def get_embedder(model_name: str, **kwargs) -> Embedder:
    """Factory function to get an embedder instance by model name and config."""
    if 'embedding-001' in model_name:
        # Import the new GeminiEmbedder
        return GeminiEmbedder(model_name=model_name, **kwargs)
    elif 'text-embedding' in model_name:
        # from .models.openai_embedder import OpenAIEmbedder
        return OpenAIEmbedder(model_name=model_name, **kwargs)
    elif model_name in HUGGINGFACE_MODELS:
        # from .models.scibert_embedder import SciBertEmbedder
        return HuggingfaceEmbedder(model_name=model_name, **kwargs)
    else: # Assume sentence-transformer
        # from .models.sentence_transformer_embedder import SentenceTransformerEmbedder
        return SentenceTransformerEmbedder(model_name=model_name, **kwargs)

def precompute_corpus_embeddings(embedder: Embedder, metadata_path: Path, cache_path: Path):
    """Generates and caches embeddings for the function metadata corpus."""
    if cache_path.exists():
        print(f"Corpus embeddings already exist at '{cache_path}'. Skipping generation.")
        return

    print(f"Corpus embeddings not found. Generating for model '{embedder.model_name}'...")
    metadata = load_json(metadata_path)
    if not metadata:
        raise ValueError(f"Could not load metadata from {metadata_path}")

    descriptions = [item.get('description', '') for item in metadata]
    
    print(f"Generating embeddings for {len(descriptions)} function descriptions...")
    embeddings = embedder.get_embeddings(descriptions)

    print(f"Saving embeddings to '{cache_path}'...")
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, 'wb') as f:
        pickle.dump(embeddings, f)
    print("Corpus embeddings generated and saved successfully.")

# --- Stage Functions (largely unchanged) ---


def run_description_generation_stage(benchmark: List[Dict], model: str, cache_path: Path, force: bool, max_workers: int = None) -> Dict[str, str]:
    """Stage 1: Generate natural language descriptions for each route in the benchmark."""
    print("\n--- Stage 1: Generating Strategy Descriptions ---")
    cached_descriptions = {} if force else load_json(cache_path)
    
    items_to_process = []
    for item in benchmark:
        route_id = item['selected_route_id']
        # if route_id not in cached_descriptions: # Only process missing items
        #      items_to_process.append(item)
    
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

def run_query_rewriting_stage(descriptions: Dict[str, str], model: str, cache_path: Path, force: bool, reasoning_budget_for_rewrite: int = 0, max_workers: int = None, ) -> Dict[str, Dict]:
    """Stage 2: Rewrite natural language descriptions into structured JSON queries."""
    print("\n--- Stage 2: Rewriting Descriptions to Queries ---")
    # cached_queries = {} if force else load_json(cache_path)
    cached_queries = load_json(cache_path)
    print(f"Loaded {len(cached_queries)} cached queries from '{cache_path}'.")
    items_to_process = []
    # for route_id, desc in descriptions.items():
        # if route_id not in cached_queries and desc:
        #     items_to_process.append((route_id, desc))
    
    if not items_to_process:
        print("All queries found in cache.")
        return cached_queries
    
    print(f"Rewriting {len(items_to_process)} new descriptions into queries...")
    
    # max_workers=5 is a common choice for LLM calls to balance concurrency and API rate limits
    with ThreadPoolExecutor(max_workers=max_workers or 5) as executor: 
        # Submit all tasks
        future_to_route = {
            executor.submit(rewrite_query, desc, model, reasoning_budget_for_rewrite): route_id 
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

def run_retrieval_stage(retriever: 'StrategyRetriever', queries: Dict[str, Dict], top_k: int, top_n_functions: Union[int, None]) -> Dict[str, List[str]]:
    print(f"  Running retrieval queries (top_n_functions={top_n_functions})...")
    retrieval_results = {}
    for route_id, query in tqdm(queries.items(), desc="Processing queries", leave=False):
        if not query or 'queries' not in query:
            retrieval_results[route_id] = []
            continue
        results = retriever.retrieve_complex(query, top_k=top_k, top_n_functions=top_n_functions)
        # The original code had a small bug here, it should be res['id'] or res['route_id']
        # Based on the context, I'll assume the key is 'route_id'
        retrieval_results[route_id] = [res['route_id'] for res in results]
    return retrieval_results

def _calculate_accuracy(results: Dict[str, List[str]], top_k_values: List[int]) -> Dict[int, float]:
    total_queries = len(results)
    if total_queries == 0:
        return {k: 0.0 for k in top_k_values}
    hits = {k: sum(1 for gt, preds in results.items() if gt in preds[:k]) for k in top_k_values}
    return {k: (hits[k] / total_queries) * 100 for k in sorted(hits.keys())}

def show_accuracy_summary(accuracies: Dict[int, float], total_queries: int, top_n_display_value: Union[int, str]):
    print(f"\n--- Evaluation Summary for (top_n={top_n_display_value}) ---")
    print(f"Total Routes Evaluated: {total_queries}")
    for k, acc in sorted(accuracies.items()):
        print(f"  Top-{k} Accuracy: {acc:.2f}%")
    print("--------------------------")


# MODIFIED function
def run_retrieval_sweep_and_save_results(
    retriever: 'StrategyRetriever',
    queries: Dict[str, Dict],
    top_k_values: List[int],
    top_n_functions_sweep_values: List[Union[int, None]],
    output_dir: Path,
    embedding_config_name: str,
):
    """
    Runs a retrieval sweep over different top_n_functions values, calculates accuracy,
    and saves both a summary CSV and detailed JSON results for each run.
    """
    if not queries:
        print("No queries available for retrieval. Skipping sweep.")
        return

    # Ensure the main output and detailed results directories exist
    output_dir.mkdir(parents=True, exist_ok=True)
    detailed_results_dir = output_dir / "results"
    detailed_results_dir.mkdir(parents=True, exist_ok=True)

    # --- Setup for summary CSV file ---
    csv_header = ["top_n_functions"] + [f"Top-{k} Accuracy" for k in sorted(top_k_values)] + ["time"]
    csv_filename = output_dir / f"retrieval_sweep_{embedding_config_name}.csv"
    
    max_k_for_retrieval = max(top_k_values) if top_k_values else 1
    
    with open(csv_filename, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=csv_header)
        writer.writeheader()

        for top_n in top_n_functions_sweep_values:
            start_time = time.time()
            
            # 1. Run the retrieval for the current top_n value
            results = run_retrieval_stage(retriever, queries, top_k=max_k_for_retrieval, top_n_functions=top_n)
            
            end_time = time.time()
            
            # 2. *** NEW: Save detailed JSON results for this specific run ***
            top_n_str = 'all' if top_n is None else str(top_n)
            json_filename = detailed_results_dir / f"retrieval_results_{embedding_config_name}_topn_{top_n_str}.json"
            
            with open(json_filename, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2)
            print(f"  Detailed results saved to '{json_filename}'")
            
            # 3. Calculate and display accuracy summary (existing logic)
            accuracies = _calculate_accuracy(results, top_k_values)
            show_accuracy_summary(accuracies, len(queries), top_n)
            
            # 4. Write the summary row to the CSV file (existing logic)
            row = {"top_n_functions": top_n, "time": end_time - start_time}
            for k, acc in accuracies.items():
                row[f"Top-{k} Accuracy"] = f"{acc:.2f}"
            writer.writerow(row)
    
    print(f"\nSweep summary saved to '{csv_filename}'")
    print(f"Detailed JSON results saved in '{detailed_results_dir}'")


async def main():
    parser = argparse.ArgumentParser(description="End-to-end strategy retrieval evaluation pipeline with embedding model sweep.")
    parser.add_argument("--benchmark_path", type=str, default="data/benchmark/eval_runs/merged_bench.json", help="Path to the benchmark JSON file.")
    parser.add_argument("--metadata_path", type=str, default="data/retrieval/function_metadata_database.json", help="Path to function metadata.")
    parser.add_argument("--route_db_dir", type=str, default="data/benchmark/uspto_trees_structured", help="Directory containing route data.")
    parser.add_argument("--output_dir", type=str, default="data/benchmark/eval_runs_uspto", help="Directory to save caches and final results.")
    parser.add_argument("--description_model", type=str, default="google/gemini-2.5-pro")
    parser.add_argument("--rewriter_model", type=str, default="google/gemini-2.5-pro")
    parser.add_argument("--top_k", nargs='+', type=int, default=[1, 3, 5, 10, 20, 50])
    parser.add_argument("--top_n_functions_sweep", nargs='+', type=int, default=[40])
    parser.add_argument("--force_description",type=bool, default=False)
    parser.add_argument("--force_rewrite",type=bool, default=False)
    parser.add_argument("--reasoning_budget_for_rewrite", type=int, default=8000)
    args = parser.parse_args()

    base_path = Path(__file__).resolve().parents[2] # Simplified for standalone execution
    benchmark_path = base_path / args.benchmark_path
    metadata_path = base_path / args.metadata_path
    route_db_dir = base_path / args.route_db_dir
    output_dir = base_path / args.output_dir
    corpus_embedding_cache_dir = base_path / "data/retrieval_embeddings_uspto"
    corpus_embedding_cache_dir.mkdir(parents=True, exist_ok=True)
    
    description_cache_path = output_dir / "merged_descriptions.json"
    rewritten_query_cache_path = output_dir / f"queries/queries_google_gemini-2.5-pro.json"

    benchmark_data = load_json(benchmark_path)
    if not benchmark_data:
        print(f"Error: Benchmark file not found or empty at '{benchmark_path}'")
        return

    # Stages 1 & 2 run only once as they are model-agnostic
    descriptions = run_description_generation_stage(benchmark_data, args.description_model, description_cache_path, args.force_description)
    queries = run_query_rewriting_stage(descriptions, args.rewriter_model, rewritten_query_cache_path, args.force_rewrite, args.reasoning_budget_for_rewrite)
    
    ablation_results_dir = output_dir / "ablations"
    ablation_results_dir.mkdir(parents=True, exist_ok=True)
    
    # --- MAIN SWEEP LOOP ---
    for config in EMBEDDING_CONFIGS:
        config_name = get_config_name(config)
        print(f"\n{'='*25} RUNNING SWEEP FOR: {config_name} {'='*25}")

        # 1. Setup Embedders (Document and Query)
        doc_embedder_config = config.copy()
        query_embedder_config = config.copy()

        # Handle asymmetric Gemini retrieval case
        if config.get("task_type") == "RETRIEVAL":
            print("  Asymmetric retrieval mode detected for Gemini.")
            doc_embedder_config["task_type"] = config["doc_task_type"]
            query_embedder_config["task_type"] = config["query_task_type"]
        
        doc_embedder = get_embedder(**doc_embedder_config)
        query_embedder = get_embedder(**query_embedder_config)

        # 2. Pre-compute Corpus Embeddings using the document embedder
        corpus_embedding_cache_path = corpus_embedding_cache_dir / f"{config_name}.pck"
        precompute_corpus_embeddings(doc_embedder, metadata_path, corpus_embedding_cache_path)

        # 3. Initialize Retriever with the query embedder
        retriever = StrategyRetriever(
            metadata_db_path=str(metadata_path),
            route_db_dir=str(route_db_dir),
            embedding_cache_path=str(corpus_embedding_cache_path),
            embedder=query_embedder, # Use the query-specific embedder
        )

        # 4. Run Retrieval Sweep and Save Results
        run_retrieval_sweep_and_save_results(
            retriever=retriever,
            queries=queries,
            top_k_values=args.top_k,
            top_n_functions_sweep_values=args.top_n_functions_sweep,
            output_dir=ablation_results_dir,
            embedding_config_name=config_name,
        )

    print(f"\nAll embedding model sweeps completed. Results saved to '{ablation_results_dir}'")

if __name__ == "__main__":
    # Ensure necessary API keys (e.g., OPENAI_API_KEY, GOOGLE_API_KEY) are set as environment variables.
    asyncio.run(main())