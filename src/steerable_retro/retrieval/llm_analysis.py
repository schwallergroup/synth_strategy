# src/steerable_retro/llm_utils.py

import json
import litellm
import networkx as nx
import importlib
import requests
import os
import ast
from aizynthfinder.chem import FixedRetroReaction
from aizynthfinder.reactiontree import ReactionTree
from typing import Dict, Any

OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY")
# Configure litellm to be quiet
litellm.suppress_debug_info = True

def get_smiles_with_depth(tree: ReactionTree):
    """Get all reaction SMILES from a tree, with their depth."""
    smiles_with_depth = []
    for node in tree.graph.nodes():
        if isinstance(node, FixedRetroReaction):
            forward_smi = node.metadata["smiles"]
            
            # Get distance of node from root and correct for molecule nodes
            depth = nx.shortest_path_length(tree.graph, source=tree.root, target=node)
            depth = int((depth - 1) / 2)
            smiles_with_depth.append((depth, forward_smi))
    # Sort by depth to present the route logically
    smiles_with_depth.sort(key=lambda x: x[0])
    return smiles_with_depth

def linearize_route(route_data: Dict[str, Any]) -> str:
    """Converts a route dictionary into a linearized string representation."""
    tree = ReactionTree.from_dict(route_data)
    smiles_list = get_smiles_with_depth(tree)
    
    linearized = []
    for i, (depth, smi) in enumerate(smiles_list):
        linearized.append(f"Step {i+1} (Depth {depth}): {smi}")
    # print(f"Linearized route:\n{linearized}")
    return "\n".join(linearized)

def _call_llm(prompt: str, model: str) -> str:
    """A helper function to call the LLM API via litellm."""
    try:
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENROUTER_API_KEY}"},
            data=json.dumps({
                "model": "google/gemini-2.5-pro",
                "messages": [{"role": "user", "content": prompt}],
                "temperature": 0.1,
                "response_format": {"type": "json_object"},
                # "reasoning": {
                #                     "max_tokens": 8000, 
                #                 }
            }),
            timeout=180
        )

        response.raise_for_status()
        api_response = response.json()
        return api_response['choices'][0]['message']['content']
    except requests.exceptions.RequestException as e:
        print(f"Error calling API: {e}")
        print(f"input prompt: {prompt}")
        print("Status code:", response.status_code)
        print("Raw response:", api_response) 
        raise Exception
        return ""
    except (KeyError, IndexError) as e:
        print(f"Error parsing API response structure: {e}")
        raise Exception
        return ""

def generate_strategy_description(route_data: Dict[str, Any], model: str) -> str:
    """Generates a natural language strategy description for a given route."""
    # linearized = linearize_route(route_data)
    inputs = importlib.import_module(
                "steerable_retro.llm.prompts.strategy_description"
            ).strategy_describer
    prompt = inputs.format(synthesis_route=route_data)
    description = _call_llm(prompt, model)
    return description.strip()

def rewrite_query(description: str, model: str) -> Dict[str, Any]:
    """Rewrites a natural language description into a structured JSON query."""
    inputs = importlib.import_module(
                "steerable_retro.llm.prompts.query_rewriting_2"
            ).query_rewriting
    #" \n ".join(ast.literal_eval(description))
    prompt = inputs.replace("{natural_language_description}", " \n ".join(description))
    response_text = _call_llm(prompt, model)
    try:
        return json.loads(response_text)
    except json.JSONDecodeError:
        print(f"Error: Failed to parse JSON from LLM response:\n{response_text}")
        raise Exception
        return {}