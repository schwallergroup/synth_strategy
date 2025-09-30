import json
import os
import requests
import re
import time
from datetime import datetime, timedelta
from typing import List, Dict, Any, Tuple, Optional

# --- NEW: Concurrency and Threading Imports ---
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import tqdm  # For the progress bar (pip install tqdm)
import ast
from .llm_utils import extract_functions

# --- OpenRouter API Configuration ---
OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY") # A better way to get the key
if not OPENROUTER_API_KEY:
    raise ValueError("OPENROUTER_API_KEY environment variable not set. Please set your OpenRouter API key.")

# --- Concurrency and Rate Limit Configuration ---
MAX_CONCURRENT_REQUESTS = 20



FUNCTION_ANNOTATION_PROMPT = """
You are an expert Python programmer and a PhD-level synthetic chemist. Your task is to analyze a Python function that encodes a specific chemical synthesis strategy and generate a structured JSON metadata object that describes its capabilities.

Read the provided Python function carefully. Your goal is to deconstruct its logic into a set of atomic checks and structural constraints.

The function you are analyzing uses an external `checker` object with the following methods:
- `checker.check_ring("ring_name", "smiles_string")`: Checks for the presence of a ring system (e.g., "pyrazole", "benzene").
- `checker.check_reaction("reaction_name", "reaction_smiles")`: Checks if a reaction SMILES matches a named reaction (e.g., "Suzuki coupling").
- `checker.check_functional_group("group_name", "smiles_string")`: Checks for a functional group (e.g., "carboxylic acid").

**Your analysis should follow these steps:**

1.  **Extract the Core Description**: Read the function's docstring and summarize its high-level purpose in the `description` field.
2.  **Identify Atomic Checks**: Analyze the code to find all specific chemical entities the function explicitly checks for.
    *   **`named_reactions`**: List all string literals passed to `checker.check_reaction()`. If the function forms a ring *de novo* (e.g., pyrazole is in product but not reactants), add a "ring_formation" event. If a ring is in reactants but not product, add a "ring_destruction" event.
    *   **`ring_systems`**: List all string literals passed to `checker.check_ring()`.
    *   **`functional_groups`**: List all string literals passed to `checker.check_functional_group()`.
3.  **Identify Structural Constraints**: Analyze the control flow of the function to understand how it uses the results of the atomic checks. This captures the *logic* of the strategy.
    *   **`positional`**: If the code checks the `depth` variable, describe the constraint. Use values like "last_stage" (for `depth == 0`), "first_stage" (for `depth == max_depth`), or "not_last_stage" (for `depth > 0`).
    *   **`count`**: If the code counts occurrences of an event (e.g., `modifications >= 2`), describe the target, the comparison operator (e.g., "==", ">=", "<"), and the value.
    *   **`sequence`**: If the function checks for one event happening before or after another (e.g., a Suzuki coupling happening after an amide formation), describe the sequence.
    *   **`co-occurrence`**: If the function checks that two or more distinct events must happen anywhere in the route (e.g., `found_suzuki and found_amide`), list the required events.
    *   **`negation`**: If the function returns True only if a specific event does *not* happen, describe the negated event.

**Output Format:**

You MUST respond with a single, valid JSON object. Do not include any text, explanations, or code blocks outside of the JSON structure.

The JSON object must have the following structure:
{
  "function_id": "string", // Placeholder, you will fill this in post-processing.
  "filepath": "string", // Placeholder, you will fill this in post-processing.
  "description": "string", // High-level summary from the docstring.
  "atomic_checks": {
    "named_reactions": ["string", ...], // List of reaction names or events like "ring_formation".
    "ring_systems": ["string", ...], // List of ring systems.
    "functional_groups": ["string", ...] // List of functional groups.
  },
  "structural_constraints": [
    {
      "type": "string", // e.g., "positional", "count", "sequence", "co-occurrence", "negation"
      "details": {} // A dictionary containing specific details for the constraint type.
    },
    ...
  ]
}

**Example of a good `structural_constraints.details` object:**
- For `count`: `{"target": "pyrazole_modification", "operator": ">=", "value": 2}`
- For `positional`: `{"target": "Friedel-Crafts acylation", "position": "last_stage"}`
- For `co-occurrence`: `{"targets": ["Suzuki coupling", "amide_formation"]}`

Now, analyze the following Python function and generate the JSON metadata object.
"""


# --- Helper Functions ---

def load_function_names_from_directory(directory_path: str) -> List[str]:
    """Scans a directory for .py files and returns their base names."""
    print(f"Scanning directory for function files: {directory_path}")
    try:
        if not os.path.isdir(directory_path):
            raise FileNotFoundError(f"The specified directory does not exist: {directory_path}")
        function_names = [os.path.splitext(f)[0] for f in os.listdir(directory_path) if f.endswith('.py')]
        print(f"Found {len(function_names)} functions to process.")
        return function_names
    except FileNotFoundError as e:
        print(f"FATAL ERROR: {e}")
        return []

def load_function_code(function_name: str, base_path: str) -> Optional[str]:
    """Loads the full Python file content for a given function name."""
    file_path = os.path.join(base_path, f"{function_name}.py")
    if not os.path.exists(file_path):
        print(f"  -> WARNING: Code file not found at {file_path}")
        return None
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return f.read()
    except Exception as e:
        print(f"  -> ERROR: Failed to read file {file_path}: {e}")
        return None

def call_llm_api_for_annotation(code: str) -> str:
    """Calls the specified LLM API with the annotation prompt."""
    if not OPENROUTER_API_KEY or OPENROUTER_API_KEY == "YOUR_API_KEY_HERE":
        print("FATAL ERROR: OPENROUTER_API_KEY is not set.")
        return ""
    try:
        user_prompt = f"{FUNCTION_ANNOTATION_PROMPT}\n\n```python\n{code}\n```"
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENROUTER_API_KEY}"},
            data=json.dumps({
                "model": "google/gemini-2.5-pro",
                "messages": [{"role": "user", "content": user_prompt}],
                "temperature": 0.1,
                "response_format": {"type": "json_object"}
            }),
            timeout=180
        )
        response.raise_for_status()
        api_response = response.json()
        return api_response['choices'][0]['message']['content']
    except requests.exceptions.RequestException as e:
        print(f"Error calling API: {e}")
        return ""
    except (KeyError, IndexError) as e:
        print(f"Error parsing API response structure: {e}")
        return ""

def parse_api_json_response(response: str) -> Optional[Dict[str, Any]]:
    """Parses the JSON string response from the API."""
    if not response: return None
    try:
        return json.loads(response)
    except json.JSONDecodeError:
        print(f"  -> Warning: Could not decode JSON from response. Content: {response[:300]}...")
        return None

# --- Main Processing Logic ---

def annotate_single_function(function_name: str, code_base_path: str) -> Tuple[str, Optional[Dict[str, Any]]]:
    """
    Processes a single function: loads code, sends to LLM for annotation, and parses the result.
    """
    print(f"Annotating '{function_name}'...")
    full_code = load_function_code(function_name, code_base_path)
    if not full_code:
        return function_name, None

    response_text = call_llm_api_for_annotation(full_code)
    if not response_text:
        return function_name, None

    annotation_data = parse_api_json_response(response_text)
    return function_name, annotation_data

def generate_function_metadata_db(code_dir: str, output_json_path: str):
    """
    Main function to process all functions, generate metadata, and save to a single JSON database file.
    """
    function_names = load_function_names_from_directory(code_dir)
    if not function_names:
        print("No functions to process. Exiting.")
        return

    all_annotations = []
    error_count = 0
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENT_REQUESTS) as executor:
        future_to_function = {
            executor.submit(annotate_single_function, name, code_dir): name for name in function_names
        }
        progress_bar = tqdm.tqdm(as_completed(future_to_function), total=len(function_names), desc="Generating Function Metadata")
        
        for future in progress_bar:
            try:
                name, result_data = future.result()
                if result_data:
                    # Enrich the LLM's response with file system info
                    result_data['function_id'] = name
                    result_data['filepath'] = os.path.join(code_dir, f"{name}.py")
                    all_annotations.append(result_data)
                else:
                    error_count += 1
                    print(f"Failed to get valid annotation for '{name}'.")
            except Exception as exc:
                name = future_to_function[future]
                print(f"'{name}' generated an exception: {exc}")
                error_count += 1

    print(f"\nAll tasks completed. Saving metadata database to {output_json_path}...")
    # Ensure the output directory exists
    os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
    with open(output_json_path, 'w') as f:
        json.dump(all_annotations, f, indent=2)
    
    print("\n--- Final Summary ---")
    print(f"Total functions processed: {len(function_names)}")
    print(f"Successfully annotated:    {len(all_annotations)}")
    print(f"Errors during annotation:  {error_count}")
    print(f"Metadata database saved to '{output_json_path}'")
    print("-----------------------------")
