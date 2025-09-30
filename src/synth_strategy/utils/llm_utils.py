# src/steerable_retro/llm_utils.py

import json
import litellm
import networkx as nx
import importlib
import requests
import os
import ast
import contextlib
import io
import re

from synth_strategy.utils import check
# import ast # Not used anymore
from aizynthfinder.chem import FixedRetroReaction
from aizynthfinder.reactiontree import ReactionTree
from typing import Dict, Any, Optional, Tuple, List # Added Optional

OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY")
# Configure litellm to be quiet
litellm.suppress_debug_info = True

REASONING_BY_DEFAULT = {"google/gemini-2.5-pro", "anthropic/claude-sonnet-3.7:thinking"}
REASONING_BY_TOGGLE_GOOGLE = {"google/gemini-2.5-flash", "anthropic/claude-sonnet-4", "anthropic/claude-sonnet-3.7"}

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

# Modified: Added optional reasoning_budget parameter
def _call_llm(prompt: str, model: str, reasoning_budget: Optional[int] = None) -> str:
    """A helper function to call the LLM API via litellm."""
    response = None # Initialize response to None for error handling
    api_response = None # Initialize api_response for error handling
    try:
        payload = {
            "model": model,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": 0.1,
            "response_format": {"type": "json_object"},
        }

        if model in REASONING_BY_DEFAULT:
            # No specific reasoning toggle or budget for these models
            pass
        elif model in REASONING_BY_TOGGLE_GOOGLE:
            # Apply reasoning budget if provided, otherwise use default of 8000 tokens
            payload["reasoning"] = {"max_tokens": reasoning_budget if reasoning_budget is not None else 0}
        
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENROUTER_API_KEY}"},
            data=json.dumps(payload),
            timeout=180
        )

        response.raise_for_status() # Raise HTTPError for bad responses (4xx or 5xx)
        api_response = response.json()
        return api_response['choices'][0]['message']['content']
    except requests.exceptions.RequestException as e:
        print(f"Error calling API: {e}")
        print(f"Input prompt: {prompt}")
        if response is not None:
            print("Status code:", response.status_code)
            print("Raw response body:", response.text)
        else:
            print("No HTTP response received.")
        raise Exception("LLM API call failed.")
    except (KeyError, IndexError) as e:
        print(f"Error parsing API response structure: {e}")
        if api_response is not None:
            print("Raw API response:", api_response)
        else:
            print("API response was not a valid JSON or missing expected keys.")
        raise Exception("LLM API response parsing failed.")

def generate_strategy_description(route_data: Dict[str, Any], model: str) -> str:
    """Generates a natural language strategy description for a given route."""
    inputs = importlib.import_module(
                "synth_strategy.llm.prompts.strategy_description"
            ).strategy_describer
    prompt = inputs.format(synthesis_route=route_data)
    # Reasoning budget is not applied to description generation phase
    description = _call_llm(prompt, model) 
    return description.strip()

# Modified: Added optional reasoning_budget_for_rewrite parameter
def rewrite_query(description: str, model: str, reasoning_budget_for_rewrite: Optional[int] = None) -> Dict[str, Any]:
    """Rewrites a natural language description into a structured JSON query."""
    inputs = importlib.import_module(
                "synth_strategy.llm.prompts.query_rewriting_2"
            ).query_rewriting
    
    # Corrected: `description` is expected to be a single string from `generate_strategy_description`
    # The original commented-out `ast.literal_eval` and `join` were likely remnants of a different approach.
    prompt = inputs.replace("{natural_language_description}", description)
    
    # Pass reasoning_budget_for_rewrite to _call_llm
    response_text = _call_llm(prompt, model, reasoning_budget=reasoning_budget_for_rewrite)
    try:
        return json.loads(response_text)
    except json.JSONDecodeError:
        print(f"Error: Failed to parse JSON from LLM response:\n{response_text}")
        raise Exception("Failed to parse JSON from LLM query rewrite response.")
    


def extract_imports(code_str: str) -> str:
    """
    Parses the code string using ast and returns a string with all
    top–level import (or from…import) statements.
    """
    lines = code_str.splitlines()
    import_lines = []
    try:
        module = ast.parse(code_str)
    except Exception as e:
        # If parsing fails return an empty string.
        return ""
    for node in module.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            # Use node.lineno and node.end_lineno (Python 3.8+) to capture the lines.
            start = node.lineno - 1
            end = getattr(node, "end_lineno", node.lineno)
            import_lines.extend(lines[start:end])
    return "\n".join(import_lines)


def extract_functions(code_str: str) -> list:
    """
    Parses the code string and returns a list of tuples (func_name, func_code)
    for each top–level function defined in the code.
    """
    funcs: list = []
    try:
        module = ast.parse(code_str)
    except Exception:
        return funcs
    lines = code_str.splitlines()
    for node in module.body:
        if isinstance(node, ast.FunctionDef):
            # node.lineno is 1-based.
            start = node.lineno - 1
            # node.end_lineno is available in Python 3.8+; fallback to grabbing just the header if not.
            end = getattr(node, "end_lineno", start + 1)
            func_code = "\n".join(lines[start:end])
            funcs.append((node.name, func_code))
    return funcs


def test_generated_code_and_capture(
    route_data: dict, code_str: str, checker
) -> Tuple[bool, str, str]:
    """
    Executes the code snippet in a restricted local namespace,
    expects a callable that takes a single argument (route_data).
    Provides access to a Check object as 'checker' in the candidate function.
    Returns: (passed, stdout, error_message)
    """
    # Create a namespace for exec, and insert the checker under 'checker'
    env = {"checker": checker}
    captured_stdout = ""
    error_message = ""
    passed = False
    # try:
    # Redirect stdout so that prints in the executed code will be captured
    with contextlib.redirect_stdout(io.StringIO()) as f:
        try:
            exec(code_str, env)
            # Attempt to get a function called "test_route"
            candidate_func = env.get("test_route")
            # If not found, search for a callable defined by the user.
            if not callable(candidate_func):
                for name, obj in env.items():
                    if callable(obj) and not name.startswith("__"):
                        candidate_func = obj
                        break
            if candidate_func and callable(candidate_func):
                # Call the candidate function using the provided route_data
                result = candidate_func(route_data)
                passed = bool(result)
            else:
                error_message = (
                    "No valid test function found in executed code."
                )
        except Exception as e:
            error_message = str(e)
            passed = False

        captured_stdout = f.getvalue()
    # except Exception as e:
    #     error_message = str(e)
    return passed, captured_stdout, error_message


def parse_generated_code(response_str: str) -> List[str]:
    """
    Find all Python code segments between <generated_code> and </generated_code>.
    Returns them as a list of strings.
    """
    # I want this to also work with code where its just in ```python ... ```
    if "<code>" not in response_str:
        # Look for ```python ... ``` blocks
        matches = re.findall(r"```python(.*?)```", response_str, flags=re.DOTALL)
        return [matches[-1].strip()] if matches else []
    matches = re.findall(r"<code>(.*?)</code>", response_str, flags=re.DOTALL)
    return matches


def parse_extractable_strategy(response_str: str) -> str:
    """
    Find the synthetic strategy within <synthetic_strategy>...</synthetic_strategy> tags.
    If not found, return an empty string.
    """
    match = re.search(
        r"<extractable_strategy>(.*?)</extractable_strategy>",
        response_str,
        flags=re.DOTALL,
    )
    if match:
        return match.group(1)
    return ""
