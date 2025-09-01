import json
import numpy as np
from openai import OpenAI
from sklearn.metrics.pairwise import cosine_similarity
from typing import List, Dict, Any
import os
import pickle

class HybridStrategyRetriever:
    """
    A hybrid retrieval system that uses semantic search on descriptions and
    categorical matching on specific chemical entities with partial match scoring.
    Now uses OpenAI's text-embedding-ada-002 model for embeddings.
    """
    def __init__(self,
                 metadata_db_path: str,
                 route_db_dir: str,
                 embedding_cache_path: str,
                 embedder: Embedder,
                 use_structural_constraints: bool = False,
                 cache_dir: str = ".retriever_cache"):
        print("Initializing StrategyRetriever...")
        
        self.embedder = embedder
        self.use_structural_constraints = use_structural_constraints
        self.route_db_dir = Path(route_db_dir)
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        db_dir_hash = str(hash(str(self.route_db_dir.resolve())))
        self.index_cache_path = self.cache_dir / f"inverted_index_{db_dir_hash}.pkl"

        print("Loading function metadata...")
        self.func_metadata = self._load_metadata(metadata_db_path)
        if not self.func_metadata:
            raise ValueError("Failed to load metadata or database is empty.")
        
        # Store original mappings temporarily for filtering purposes
        original_func_id_to_metadata = {func['function_id'] + ".py": func for func in self.func_metadata}
        
        print(f"Loading pre-computed description embeddings from '{embedding_cache_path}'...")
        self.embedding_cache_path = Path(embedding_cache_path)
        if not self.embedding_cache_path.exists():
            raise FileNotFoundError(
                f"Embedding cache file not found at '{self.embedding_cache_path}'. "
                "Please generate it first using the pre-computation step in the evaluation pipeline."
            )
        with open(self.embedding_cache_path, "rb") as f:
            self.description_embeddings = np.array(pickle.load(f))
        
        # Store original func_ids and its index mapping temporarily
        original_func_ids = list(original_func_id_to_metadata.keys())
        original_func_id_to_embedding_idx = {fid: i for i, fid in enumerate(original_func_ids)}

        print("Creating sets for fast atomic lookups (pre-filtering)...")
        # Store original metadata_atomic_sets temporarily
        original_metadata_atomic_sets = self._create_atomic_check_sets(self.func_metadata)
        
        # Load route database and create index. This must happen before filtering,
        # as we need the index to determine which functions have "passed".
        self.route_database, self.function_to_routes_index = self._load_or_create_index()
        if not self.route_database:
            raise ValueError("Failed to load route data or directory is empty.")
        print(f"Loaded {len(self.route_database)} routes and {len(self.function_to_routes_index)} functions in index.")

        # --- NEW LOGIC START: Filter all data structures based on passing functions ---
        print("Filtering functions to include only those passing at least once in the route database...")

        # 1. Identify all unique function IDs that are part of any 'passing_functions'
        # The function_to_routes_index already contains only these.
        all_passing_func_ids = set(self.function_to_routes_index.keys())

        if not all_passing_func_ids:
            raise ValueError(
                "No functions found that pass at least once in the route database. "
                "This might indicate an issue with route data consistency, "
                "or your metadata/embedding files contain no passing functions."
            )

        # Temporary lists to build the filtered data
        new_func_ids = []
        new_description_embeddings_list = []
        new_func_metadata_list = []
        new_metadata_atomic_sets_list = []

        # Iterate through the original `self.func_metadata`, which is the primary source
        # This ensures proper alignment between metadata, func_ids, and atomic sets.
        for i, func_meta in enumerate(self.func_metadata):
            func_id = func_meta['function_id'] + ".py"
            if func_id in all_passing_func_ids:
                new_func_ids.append(func_id)
                
                # Retrieve embedding using its original index
                original_emb_idx = original_func_id_to_embedding_idx.get(func_id)
                if original_emb_idx is not None and 0 <= original_emb_idx < len(self.description_embeddings):
                    new_description_embeddings_list.append(self.description_embeddings[original_emb_idx])
                else:
                    # Fallback for robustness, though this case should ideally not happen
                    print(f"Warning: Embedding for func_id '{func_id}' not found or out of bounds. Using a zero vector.")
                    # Ensure a zero vector of correct dimension is added
                    new_description_embeddings_list.append(np.zeros(self.description_embeddings.shape[1]))

                new_func_metadata_list.append(func_meta)
                
                # Retrieve atomic sets using its original index
                if 0 <= i < len(original_metadata_atomic_sets):
                    new_metadata_atomic_sets_list.append(original_metadata_atomic_sets[i])
                else:
                    # Fallback for robustness
                    print(f"Warning: Atomic sets for func_id '{func_id}' not found at original index {i}. Using empty dict.")
                    new_metadata_atomic_sets_list.append({})

        # 4. Update the instance attributes with the filtered data
        self.func_ids = new_func_ids
        self.description_embeddings = np.array(new_description_embeddings_list)
        self.func_metadata = new_func_metadata_list
        self.metadata_atomic_sets = new_metadata_atomic_sets_list

        # 5. Rebuild all mappings based on the new, filtered, and re-indexed lists
        self.func_id_to_embedding_idx = {fid: i for i, fid in enumerate(self.func_ids)}
        self.func_id_to_metadata = {func['function_id'] + ".py": func for func in self.func_metadata}
        self.func_id_to_atomic_sets = {
            self.func_ids[i]: self.metadata_atomic_sets[i]
            for i in range(len(self.func_ids))
        }
        
        print(f"Filtered to {len(self.func_ids)} functions that passed at least once and updated all internal structures.")
        # --- NEW LOGIC END ---

        print("Retriever initialized successfully.")

    def _load_metadata(self, db_path: str) -> List[Dict[str, Any]]:
        try:
            with open(db_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            print(f"Error loading metadata database: {e}")
            return []

    def _get_embeddings(self, texts: List[str], batch_size: int = 100) -> np.ndarray:
        """
        Get embeddings from OpenAI's text-embedding-ada-002 model.
        Processes in batches to handle API rate limits.
        
        Args:
            texts (List[str]): List of texts to embed
            batch_size (int): Number of texts to process in each batch
            
        Returns:
            np.ndarray: Array of embeddings
        """
        embeddings = []
        
        for i in range(0, len(texts), batch_size):
            batch = texts[i:i + batch_size]
            print(f"Processing embedding batch {i//batch_size + 1}/{(len(texts)-1)//batch_size + 1}")
            
            try:
                response = self.client.embeddings.create(
                    input=batch,
                    model="text-embedding-ada-002"
                )
                
                batch_embeddings = [item.embedding for item in response.data]
                embeddings.extend(batch_embeddings)
                
            except Exception as e:
                print(f"Error getting embeddings for batch starting at index {i}: {e}")
                # Add zero embeddings as fallback
                embeddings.extend([[0.0] * 1536] * len(batch))  # ada-002 has 1536 dimensions
        
        return np.array(embeddings)

    def _get_single_embedding(self, text: str) -> np.ndarray:
        """
        Get embedding for a single text using OpenAI's text-embedding-ada-002.
        
        Args:
            text (str): Text to embed
            
        Returns:
            np.ndarray: Embedding vector
        """
        try:
            response = self.client.embeddings.create(
                input=text,
                model="text-embedding-ada-002"
            )
            return np.array(response.data[0].embedding)
        except Exception as e:
            print(f"Error getting embedding for text: {e}")
            return np.zeros(1536)  # Return zero vector as fallback

    def _create_atomic_check_sets(self) -> List[Dict[str, set]]:
        """Converts atomic check lists to sets for efficient lookups."""
        check_sets = []
        for func in self.metadata:
            checks = func.get('atomic_checks', {})
            check_sets.append({
                'named_reactions': set(item.lower() for item in checks.get('named_reactions', [])),
                'ring_systems': set(item.lower() for item in checks.get('ring_systems', [])),
                'functional_groups': set(item.lower() for item in checks.get('functional_groups', []))
            })
        return check_sets

    def _get_categorical_score(self, func_idx: int, required_entities: Dict[str, List[str]]) -> float:
        """
        Calculates the categorical match score using partial matching.
        
        Scoring system:
        - Each matched entity contributes to the score
        - Score is normalized by the total number of required entities
        - Returns a value between 0.0 and 1.0
        
        Args:
            func_idx (int): Index of the function to check
            required_entities (Dict[str, List[str]]): Required entities by category
            
        Returns:
            float: Categorical score between 0.0 and 1.0
        """
        func_sets = self.atomic_check_sets[func_idx]
        
        total_required = 0
        total_matches = 0
        
        for category, required_list in required_entities.items():
            if not required_list:
                continue
                
            # Count total required entities
            total_required += len(required_list)
            
            # Normalize required items to lower case
            required_set = set(item.lower() for item in required_list)
            
            # Count matches in this category
            if func_sets.get(category):
                matches = len(required_set.intersection(func_sets[category]))
                total_matches += matches
        
        # Avoid division by zero
        if total_required == 0:
            return 0.0
            
        # Return normalized score (0.0 to 1.0)
        return total_matches / total_required

    def _get_categorical_score_with_weights(self, func_idx: int, required_entities: Dict[str, List[str]], 
                                          category_weights: Dict[str, float] = None) -> float:
        """
        Alternative categorical scoring with category weights.
        
        Args:
            func_idx (int): Index of the function to check
            required_entities (Dict[str, List[str]]): Required entities by category
            category_weights (Dict[str, float], optional): Weights for different categories
            
        Returns:
            float: Weighted categorical score between 0.0 and 1.0
        """
        if category_weights is None:
            category_weights = {
                'named_reactions': 1.0,
                'ring_systems': 0.8,
                'functional_groups': 0.6
            }
        
        func_sets = self.atomic_check_sets[func_idx]
        
        total_weighted_score = 0.0
        total_possible_score = 0.0
        
        for category, required_list in required_entities.items():
            if not required_list:
                continue
                
            weight = category_weights.get(category, 1.0)
            
            # Normalize required items to lower case
            required_set = set(item.lower() for item in required_list)
            
            # Calculate category score
            if func_sets.get(category):
                matches = len(required_set.intersection(func_sets[category]))
                category_score = matches / len(required_list)  # 0.0 to 1.0
            else:
                category_score = 0.0
            
            total_weighted_score += category_score * weight
            total_possible_score += weight
        
        # Avoid division by zero
        if total_possible_score == 0:
            return 0.0
            
        return total_weighted_score / total_possible_score

    def retrieve(self, query_strategy: Dict[str, Any], top_k: int = 3, 
                use_weighted_categorical: bool = False, 
                category_weights: Dict[str, float] = None) -> List[Dict[str, Any]]:
        """
        Retrieves the top_k most relevant functions for a single structured strategy query.
        
        Args:
            query_strategy (Dict[str, Any]): The query strategy
            top_k (int): Number of results to return
            use_weighted_categorical (bool): Whether to use weighted categorical scoring
            category_weights (Dict[str, float], optional): Category weights for weighted scoring
        
        Returns:
            List[Dict[str, Any]]: Retrieved functions with scores
        """
        query_description = query_strategy.get("description", "")
        required_entities = query_strategy.get("required_entities", {})
        
        # --- Stage 1: Semantic Score ---
        query_embedding = self._get_single_embedding(query_description)
        
        # Compute cosine similarity using sklearn
        semantic_scores = cosine_similarity([query_embedding], self.description_embeddings)[0]

        # --- Stage 2: Enhanced Categorical Scoring ---
        results = []
        for i, func in enumerate(self.metadata):
            semantic_score = semantic_scores[i]
            
            # Skip functions with very low semantic similarity
            if semantic_score < 0.1:  # Lowered threshold for partial matches
                continue
            
            # Choose categorical scoring method
            if use_weighted_categorical:
                categorical_score = self._get_categorical_score_with_weights(i, required_entities, category_weights)
            else:
                categorical_score = self._get_categorical_score(i, required_entities)
            
            # --- Hybrid Scoring Strategy ---
            # Option 1: Multiplicative (requires both semantic and categorical match)
            # final_score = semantic_score * categorical_score
            
            # Option 2: Additive with weights (allows good semantic OR good categorical)
            alpha = 0.6  # Weight for semantic score
            beta = 0.4   # Weight for categorical score
            final_score = alpha * semantic_score + beta * categorical_score
            
            # Only include results with some relevance
            if final_score > 0.1:  # Lowered threshold
                results.append({
                    'function_id': func['function_id'],
                    'final_score': float(final_score),
                    'semantic_score': float(semantic_score),
                    'categorical_score': float(categorical_score),
                    'categorical_match': categorical_score > 0,
                    'description': func['description']
                })
        
        results.sort(key=lambda x: x['final_score'], reverse=True)
        return results[:top_k]

    def retrieve_debug(self, query_strategy: Dict[str, Any], top_k: int = 5) -> List[Dict[str, Any]]:
        """
        Debug version that shows detailed scoring information.
        """
        query_description = query_strategy.get("description", "")
        required_entities = query_strategy.get("required_entities", {})
        
        print(f"\nDEBUG: Query - {query_description}")
        print(f"DEBUG: Required entities - {required_entities}")
        
        query_embedding = self._get_single_embedding(query_description)
        semantic_scores = cosine_similarity([query_embedding], self.description_embeddings)[0]

        results = []
        for i, func in enumerate(self.metadata):
            semantic_score = semantic_scores[i]
            categorical_score = self._get_categorical_score(i, required_entities)
            
            # Show detailed scoring for top candidates
            if semantic_score > 0.2 or categorical_score > 0.2:
                print(f"\nDEBUG: Function {func['function_id']}")
                print(f"  Semantic: {semantic_score:.3f}")
                print(f"  Categorical: {categorical_score:.3f}")
                print(f"  Available entities: {self.atomic_check_sets[i]}")
            
            final_score = 0.6 * semantic_score + 0.4 * categorical_score
            
            if final_score > 0.1:
                results.append({
                    'function_id': func['function_id'],
                    'final_score': float(final_score),
                    'semantic_score': float(semantic_score),
                    'categorical_score': float(categorical_score),
                    'description': func['description']
                })
        
        results.sort(key=lambda x: x['final_score'], reverse=True)
        return results[:top_k]

# --- Example Usage ---

# Initialize the retriever (make sure to set your OpenAI API key)
# You can either pass it directly or set the OPENAI_API_KEY environment variable
retriever = HybridStrategyRetriever(
    '../data/retrieval/function_metadata_database.json',
    openai_api_key="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"  # Or set OPENAI_API_KEY environment variable
)

# Define test queries to show partial matching
user_queries = [
    {
        "strategy_id": "early heterocycle formation",
        "description": "preserve a trifluoromethyl compound during the synthesis",
        "required_entities": {
            "named_reactions": [],
        }
    },
    # {
    #     "strategy_id": "late heterocycle into",
    #     "description": "late stage heterocycle introduction",
    #     "required_entities": {
    #         "named_reactions": [],
    #     }
    # },
    # {
    #     "strategy_id": "preserve heterocycle ",
    #     "description": "preservation of a heterocycle (tetrazole) in all stages of the synthesis from the starting material",
    #     "required_entities": {
    #         "ring_structures": [],
    #     }
    # },
    # {
    #     "strategy_id": "diversity heterocycle ",
    #     "description": "a high diversity of 3 or more heterocycles",
    #     "required_entities": {
    #         "named_reactions": [],
    #     }
    # },
]

# Test both queries to show difference in scoring
for query in user_queries:
    print(f"\n{'='*60}")
    print(f"Results for: {query['strategy_id']}")
    print(f"Description: {query['description']}")
    print(f"Required entities: {query['required_entities']}")
    print(f"{'='*60}")
    
    # Standard partial matching
    print("\n--- Standard Partial Matching ---")
    results = retriever.retrieve(query, top_k=10)


    for result in results:
        print(f"Function: {result['function_id']}")
        print(f"  Final Score: {result['final_score']:.3f}")
        print(f"  Semantic: {result['semantic_score']:.3f}")
        print(f"  Categorical: {result['categorical_score']:.3f}")
        print(f"  Description: {result['description'][:100]}...")
        print()
    
    # # Weighted partial matching
    # print("\n--- Weighted Partial Matching ---")
    # weighted_results = retriever.retrieve(query, top_k=3, use_weighted_categorical=True)
    # for result in weighted_results:
    #     print(f"Function: {result['function_id']}")
    #     print(f"  Final Score: {result['final_score']:.3f}")
    #     print(f"  Semantic: {result['semantic_score']:.3f}")
    #     print(f"  Categorical: {result['categorical_score']:.3f}")
    #     print(f"  Description: {result['description'][:100]}...")
    #     print()