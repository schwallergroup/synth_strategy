# src/steerable_retro/llm_utils.py

import json
import litellm
import networkx as nx
from aizynthfinder.chem import FixedRetroReaction
from aizynthfinder.reactiontree import ReactionTree
from typing import Dict, Any

from llm.prompts import STRATEGY_DESCRIPTION_PROMPT_TEMPLATE, QUERY_REWRITING_PROMPT_TEMPLATE

# Configure litellm to be quiet
litellm.suppress_debug_info = True

def get_smiles_with_depth(tree: ReactionTree):
    """Get all reaction SMILES from a tree, with their depth."""
    smiles_with_depth = []
    for node in tree.graph.nodes():
        if isinstance(node, FixedRetroReaction):
            rsmi = node.metadata["rsmi"].split(">")
            # Format as reactant>>product for forward synthesis direction
            forward_smi = f"{rsmi[-1]}>>{rsmi[0]}"
            
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
    return "\n".join(linearized)

async def _call_llm(prompt: str, model: str, max_tokens: int = 2048) -> str:
    """A helper function to call the LLM API via litellm."""
    try:
        response = await litellm.acompletion(
            model=model,
            messages=[{"role": "user", "content": prompt}],
            temperature=0.0,
            max_tokens=max_tokens,
        )
        return response.choices[0].message.content
    except Exception as e:
        print(f"Error calling LLM '{model}': {e}")
        return ""

async def generate_strategy_description(route_data: Dict[str, Any], model: str) -> str:
    """Generates a natural language strategy description for a given route."""
    linearized = linearize_route(route_data)
    prompt = STRATEGY_DESCRIPTION_PROMPT_TEMPLATE.format(linearized_route=linearized)
    description = await _call_llm(prompt, model)
    return description.strip()

async def rewrite_query(description: str, model: str) -> Dict[str, Any]:
    """Rewrites a natural language description into a structured JSON query."""
    prompt = QUERY_REWRITING_PROMPT_TEMPLATE.format(user_natural_language_input=description)
    response_text = await _call_llm(prompt, model)
    
    try:
        # Clean the response to extract only the JSON part
        json_str = response_text[response_text.find('{'):response_text.rfind('}')+1]
        return json.loads(json_str)
    except json.JSONDecodeError:
        print(f"Error: Failed to parse JSON from LLM response:\n{response_text}")
        return {}