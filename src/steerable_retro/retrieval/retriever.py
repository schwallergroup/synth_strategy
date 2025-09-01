# src/steerable_retro/retrieval.py

import json
import pickle
import os
import copy
from datetime import datetime
from pathlib import Path
from collections import defaultdict
from typing import Dict, List, Any, Set, Tuple

import orjson
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from steerable_retro.models import Embedder


class StrategyRetriever:
    """
    V8: Model-agnostic semantic-first retriever.
    Accepts an embedder object to allow for flexible embedding models.
    Prioritizes finding the best partial matches by creating a broad candidate pool
    and then scoring based on match count and semantic similarity.
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

    def _get_single_embedding(self, text: str) -> np.ndarray:
        """Uses the injected embedder to get a single embedding."""
        if not text:
            # Return a zero vector of the correct dimension
            return np.zeros(self.description_embeddings.shape[1])
        try:
            return self.embedder.get_embeddings([text])[0]
        except Exception as e:
            print(f"Error getting embedding for text: '{text[:50]}...'. Error: {e}")
            return np.zeros(self.description_embeddings.shape[1])

    # --- All other methods from the previous version of StrategyRetriever ---
    # --- (_find_top_n_similar_functions, _load_or_create_index, etc.)   ---
    # --- are included here without change. For brevity, they are omitted ---
    # --- but they should be copied directly from the previous version.   ---
    # --- The only change was in __init__ and _get_single_embedding.      ---

    #<-- PASTE ALL OTHER METHODS FROM THE PREVIOUS StrategyRetriever HERE -->
    def _find_top_n_similar_functions(self, nl_description: str, top_n: int) -> Set[str]:
        if not nl_description:
            return set(self.func_ids)
        query_embedding = self._get_single_embedding(nl_description)
        similarities = cosine_similarity([query_embedding], self.description_embeddings)[0]
        top_n_indices = np.argsort(similarities)[-top_n:]
        return {self.func_ids[i] for i in top_n_indices}

    # --- NEW HELPER METHOD ---
    def _find_top_n_similar_functions_with_embedding(self, query_embedding: np.ndarray, top_n: int) -> Set[str]:
        """Finds top N similar functions using a pre-computed embedding."""
        if np.all(query_embedding == 0):
            # Handles cases where the description was empty
            return set(self.func_ids)
        similarities = cosine_similarity([query_embedding], self.description_embeddings)[0]
        top_n_indices = np.argsort(similarities)[-top_n:]
        return {self.func_ids[i] for i in top_n_indices}

    def _get_latest_data_mtime(self) -> float:
        if not self.route_db_dir.is_dir(): return 0.0
        files = list(self.route_db_dir.rglob('*.json'))
        if not files: return 0.0
        return max(p.stat().st_mtime for p in files)

    def _load_or_create_index(self) -> Tuple[Dict, Dict]:
        cache_is_valid = False
        latest_data_mtime = self._get_latest_data_mtime()
        if self.index_cache_path.exists() and latest_data_mtime > 0:
            cache_mtime = self.index_cache_path.stat().st_mtime
            if cache_mtime > latest_data_mtime:
                cache_is_valid = True
        if cache_is_valid:
            print(f"Loading route database and index from cache: '{self.index_cache_path}'")
            try:
                with open(self.index_cache_path, 'rb') as f:
                    cached_data = pickle.load(f)
                return cached_data['route_db'], cached_data['index']
            except (pickle.UnpicklingError, KeyError, IOError) as e:
                print(f"Warning: Could not load cache file ({e}). Rebuilding.")
        print("Cache not found or is stale. Building from source...")
        route_db = self._load_route_data(str(self.route_db_dir))
        if not route_db: raise ValueError("Failed to load route data or directory is empty.")
        index = self._create_inverted_index(route_db)
        print(f"Saving new cache to '{self.index_cache_path}'")
        try:
            with open(self.index_cache_path, 'wb') as f:
                pickle.dump({'route_db': route_db, 'index': index}, f)
        except (IOError, pickle.PicklingError) as e:
            print(f"Warning: Could not save cache file: {e}")
        return route_db, index

    def _create_inverted_index(self, route_db: Dict[str, Dict]) -> Dict[str, Set[str]]:
        print("Creating inverted index from functions to routes...")
        index = defaultdict(set)
        for route_id, route_data in route_db.items():
            for func_id in route_data.get('passing_functions', {}).keys():
                index[func_id].add(route_id)
        print(f"Inverted index created with {len(index)} unique functions.")
        return index
    
    def merge_dicts(self, dict_list):
        # Get all unique keys
        all_keys = set().union(*(d.keys() for d in dict_list))
        
        # Combine values for each key
        result = {}
        for key in all_keys:
            values = []
            for d in dict_list:
                if key in d:
                    if isinstance(d[key], list):
                        values.extend(d[key])
                    else:
                        values.append(d[key])
            result[key] = values
        
        return result

    def _load_route_data(self, dir_path: str) -> Dict[str, Dict[str, Any]]:
        route_db = {}
        p = Path(dir_path)
        if not p.is_dir():
            print(f"Error: Route directory not found at '{dir_path}'")
            return {}
        json_files = list(p.glob('*.json'))
        print(f"Found {len(json_files)} route files to load.")
        for file_path in json_files:
            try:
                with open(file_path, 'rb') as f:
                    route_data = orjson.loads(f.read())
                    for idx, data in enumerate(route_data):
                        if isinstance(data, dict) and 'passing_functions' in data:
                            route_db[f"{file_path.stem}_{idx}"] = data
            except (orjson.JSONDecodeError, IOError) as e:
                print(f"Warning: Could not load or parse {file_path}: {e}")
        return route_db

    def _load_metadata(self, db_path: str) -> List[Dict[str, Any]]:
        try:
            with open(db_path, 'r', encoding='utf-8') as f: return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            print(f"Error loading metadata database: {e}")
            return []

    def _create_atomic_check_sets(self, metadata_list: List[Dict]) -> List[Dict[str, Set[str]]]:
        return [{
            cat: set(item.lower() for item in func_meta.get('atomic_checks', {}).get(cat, []))
            for cat in ['named_reactions', 'ring_systems', 'functional_groups']
        } for func_meta in metadata_list]

    def _evaluate_atomic_checks(self, filters: Dict, atomic_sets: Dict[str, Set[str]]) -> bool:
        # CORRECT: Handle implicit AND for a list of filters at the start.
        # A query like [{"A":...}, {"B":...}] means A AND B.
        if isinstance(filters, list):
            return all(self._evaluate_atomic_checks(sub_filter, atomic_sets) for sub_filter in filters)

        # The rest of the logic now correctly assumes 'filters' is a dictionary
        if 'AND' in filters:
            and_clause = filters['AND']
            if isinstance(and_clause, list):
                return all(self._evaluate_atomic_checks(sub_filter, atomic_sets) for sub_filter in and_clause)
            elif isinstance(and_clause, dict):
                return self._evaluate_atomic_checks(and_clause, atomic_sets)
            else:
                 raise ValueError(f"The value for an 'AND' key must be a list or a dict, not {type(and_clause)}")
        
        if 'OR' in filters:
            or_clause = filters['OR']
            if isinstance(or_clause, list):
                return any(self._evaluate_atomic_checks(sub_filter, atomic_sets) for sub_filter in or_clause)
            elif isinstance(or_clause, dict):
                for category, or_terms in or_clause.items():
                    if not isinstance(or_terms, list):
                        raise ValueError(f"In a compact OR clause, the value for category '{category}' must be a list of terms.")
                    if not {term.lower() for term in or_terms}.isdisjoint(atomic_sets.get(category, set())):
                        return True # An item in the OR clause matched, so the whole OR is true.
                return False # No items in the compact OR clause matched.
            else:
                raise ValueError(f"The value for an 'OR' key must be a list or a dict, not {type(or_clause)}")
        
        if 'NOT' in filters:
            return not self._evaluate_atomic_checks(filters['NOT'], atomic_sets)
        
        # Base Case: This is an atomic check dictionary like {"ring_systems": ["piperazine"]}
        # The faulty merge_dicts call and list check are removed from here.
        for category, required_terms in filters.items():
            if not isinstance(required_terms, list):
                return False
                #  raise ValueError(f"The value for category '{category}' must be a list of terms.")
            
            # This is where the original error occurred. It will now only receive proper term lists.
            if not {term.lower() for term in required_terms}.issubset(atomic_sets.get(category, set())):
                return False
        return True

    def _filter_candidate_functions(self, filters: Dict, initial_candidates: Set[str] = None) -> Set[str]:
        if initial_candidates is None:
            initial_candidates = set(self.func_id_to_metadata.keys())
        if not filters:
            return initial_candidates
        return {
            func_id for func_id in initial_candidates
            if self.func_id_to_atomic_sets.get(func_id) and self._evaluate_atomic_checks(filters, self.func_id_to_atomic_sets[func_id])
        }

    def _check_instance_match(self, query_filters: Dict, pass_details: Dict) -> bool:
        instance_atomic_sets = {
            cat: set(item.lower() for item in val) 
            for cat, val in pass_details.get('atomic_checks', {}).items()
        }
        return self._evaluate_atomic_checks(query_filters, instance_atomic_sets)

    def _check_constraints_match(self, query_constraints: List[Dict], pass_details: Dict) -> bool:
        if not query_constraints: return True
        instance_constraints = pass_details.get('structural_constraints', [])
        for q_const in query_constraints:
            q_type, q_details = q_const.get('type'), q_const.get('details', {})
            for i_const in instance_constraints:
                if i_const.get('type') == q_type:
                    i_details = i_const.get('details', {})
                    if q_type == 'positional':
                        q_target, q_pos = q_details.get('target', '').lower(), q_details.get('position', '').lower()
                        i_pos, i_targets = i_details.get('position', '').lower(), {t.lower() for t in i_details.get('targets', [])}
                        if q_pos == i_pos and q_target in i_targets:
                            break
            else:
                return False
        return True

    def _route_has_matching_instance(self, route_data: Dict, sub_query: Dict, candidate_function_ids: Set[str]) -> Tuple[bool, List[str]]:
        filters = sub_query.get('filters', {})
        passing_functions = route_data.get('passing_functions', {})
        relevant_func_ids = candidate_function_ids.intersection(passing_functions.keys())
        
        matching_func_ids = []
        for func_id in relevant_func_ids:
            pass_details = passing_functions[func_id]
            if not isinstance(pass_details, dict): continue
            
            if not self._check_instance_match(filters, pass_details): continue
            
            constraints_match = True
            if self.use_structural_constraints:
                constraints_query = sub_query.get('constraints', [])
                constraints_match = self._check_constraints_match(constraints_query, pass_details)
            
            if constraints_match:
                matching_func_ids.append(func_id)
        return bool(matching_func_ids), matching_func_ids
    
    def _save_retrieval_results(self, query: Dict[str, Any], results: List[Dict[str, Any]]):
        output_dir = "retrieval_tests"
        try:
            os.makedirs(output_dir, exist_ok=True)
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            filename = f"retrieval_results_{timestamp}.json"
            filepath = os.path.join(output_dir, filename)
            filtered_routes_to_save = []
            for result_item in results:
                route_id = result_item['route_id']
                original_route = self.route_database.get(route_id)
                if not original_route: continue
                filtered_route = copy.deepcopy(original_route)
                useful_func_ids = {
                    func_id
                    for match_details in result_item.get('positive_matches', {}).values()
                    for func_id in match_details.get('matching_functions', [])
                }
                filtered_route['passing_functions'] = {
                    func_id: details
                    for func_id, details in original_route.get('passing_functions', {}).items()
                    if func_id in useful_func_ids
                }
                filtered_route['retrieval_info'] = {
                    'route_id': route_id,
                    'rank_score': result_item['rank_score'],
                    'match_count': result_item['match_count'],
                    'total_sub_queries': result_item['total_sub_queries'],
                    'is_full_match': result_item['is_full_match']
                }
                filtered_routes_to_save.append(filtered_route)
            output_data = {
                "query": query,
                "total_matches": len(results),
                "matching_routes": filtered_routes_to_save
            }
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(output_data, f, indent=4, ensure_ascii=False)
            print(f"Saved {len(results)} matching routes (with filtered functions) to '{filepath}'")
        except (IOError, TypeError, OSError) as e:
            print(f"Warning: Could not save retrieval results to file: {e}")

    def retrieve_complex(self, complex_query: Dict[str, Any], top_k: int = 5, top_n_functions: int = None, debug: bool = False) -> List[Dict[str, Any]]:
        operator = complex_query.get("operator", "AND").upper()
        sub_queries = complex_query.get("queries", [])
        if not sub_queries:
            return []

        if debug:
            print(f"Executing complex query with {len(sub_queries)} sub-queries (operator: '{operator}', top_n_functions: {top_n_functions}).")

        # --- MODIFICATION 1: Pre-compute all query embeddings once ---
        print("Pre-computing embeddings for all sub-queries...") # Added for clarity
        sub_query_embeddings = {
            i: self._get_single_embedding(sq_item['query'].get('natural_language_description', ''))
            for i, sq_item in enumerate(sub_queries)
        }
        print(f" -> Done. {len(sub_query_embeddings)} embeddings computed.")

        sub_query_info = []
        for i, sq_item in enumerate(sub_queries):
            sub_query = sq_item['query']
            atomically_filtered_funcs = self._filter_candidate_functions(sub_query.get('filters', {}))
            candidate_funcs = atomically_filtered_funcs

            nl_desc = sub_query.get('natural_language_description')
            if top_n_functions and nl_desc:
                # --- MODIFICATION 2: Use the pre-computed embedding ---
                query_embedding = sub_query_embeddings[i]
                semantically_filtered_funcs = self._find_top_n_similar_functions_with_embedding(query_embedding, top_n_functions)
                candidate_funcs = atomically_filtered_funcs.intersection(semantically_filtered_funcs)

            potential_routes = set()
            for func_id in candidate_funcs:
                potential_routes.update(self.function_to_routes_index.get(func_id, set()))

            sub_query_info.append({
                "routes": potential_routes,
                "negate": sq_item.get('negate', False),
                "candidate_funcs": candidate_funcs
            })

        # ... (rest of the candidate selection logic remains the same) ...
        final_candidate_route_ids = set()
        positive_route_sets = [s['routes'] for s in sub_query_info if not s['negate']]

        if positive_route_sets:
            final_candidate_route_ids = set.union(*positive_route_sets)
        else:
            final_candidate_route_ids = set(self.route_database.keys())

        if not final_candidate_route_ids:
            return []

        processed_routes = []
        for route_id in final_candidate_route_ids:
            route_data = self.route_database[route_id]
            positive_matching_funcs = {}
            sub_query_results = []
            
            for i, sq_item in enumerate(sub_queries):
                candidate_functions_for_sq = sub_query_info[i]['candidate_funcs']
                has_match, matching_funcs = self._route_has_matching_instance(
                    route_data, sq_item['query'], candidate_functions_for_sq
                )
                
                if has_match and not sq_item.get('negate', False):
                    positive_matching_funcs[i] = matching_funcs
                
                final_result = not has_match if sq_item.get('negate', False) else has_match
                sub_query_results.append(final_result)
            
            match_count = sum(1 for r in sub_query_results if r)
            is_full_match = (operator == "AND" and all(sub_query_results)) or \
                            (operator == "OR" and any(sub_query_results))
            
            if match_count > 0:
                processed_routes.append({
                    "route_id": route_id,
                    "positive_matches": positive_matching_funcs,
                    "match_count": match_count,
                    "is_full_match": is_full_match
                })
                
        if not processed_routes:
            return []
            
        scored_results = []
        for route_info in processed_routes:
            route_id = route_info['route_id']
            positive_matches = route_info['positive_matches']
            
            total_score = 0
            num_scores = 0
            for sub_query_idx, matching_func_ids in positive_matches.items():
                nl_query = sub_queries[sub_query_idx]['query'].get('natural_language_description', '')
                if not nl_query: continue
                
                # --- MODIFICATION 3: Use the pre-computed embedding from the cache ---
                # OLD LINE: query_embedding = self._get_single_embedding(nl_query)
                query_embedding = sub_query_embeddings[sub_query_idx]
                
                sub_query_score = 0
                for func_id in matching_func_ids:
                    embedding_idx = self.func_id_to_embedding_idx.get(func_id)
                    if embedding_idx is None: continue
                    
                    func_embedding = self.description_embeddings[embedding_idx]
                    score = cosine_similarity([query_embedding], [func_embedding])[0][0]
                    sub_query_score += score
                    
                if matching_func_ids:
                    total_score += sub_query_score / len(matching_func_ids)
                    num_scores += 1
                    
            final_score = (total_score / num_scores) if num_scores > 0 else 0.0
            
            scored_results.append({
                'route_id': route_id,
                'rank_score': final_score,
                'match_count': route_info['match_count'],
                'total_sub_queries': len(sub_queries),
                'is_full_match': route_info['is_full_match'],
                'positive_matches': {
                    f"sub_query_{idx}": {
                        "matching_functions": fids,
                        "query": sub_queries[idx]['query']
                    } for idx, fids in positive_matches.items()
                }
            })
        
        scored_results.sort(key=lambda x: (x['match_count'], x['rank_score']), reverse=True)
        
        final_results = scored_results[:top_k]
            
        return final_results