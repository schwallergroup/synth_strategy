# ### File: src/synth_strategy/llm/function_filtering.py (MODIFIED) ###

import json
import os
import requests
import re
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import tqdm
from typing import List, Dict, Any, Tuple

from .llm_utils import extract_functions

OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY")
if not OPENROUTER_API_KEY:
    raise ValueError("OPENROUTER_API_KEY environment variable not set.")

MAX_CONCURRENT_REQUESTS = 50

# --- Prompts remain the same ---
GEMINI_FLASH_FILTERING_PROMPT = """
You are an expert computational chemist and a senior software engineer with deep knowledge of synthetic organic chemistry. Your task is to analyze and evaluate Python functions that have been automatically generated to identify specific synthetic strategies from chemical reaction data.

**Your evaluation process allows for three types of improvements, which can occur independently or concurrently:**

1.  **Fixing Code by Removal:** If the code is buggy, inefficient, or redundant, you can improve it by **REMOVING** any number of flawed or unnecessary lines. Redundant functional group checks are very common.
2.  **Fixing Code by Conditional Modification:** You can **MODIFY a single conditional statement** (e.g., `if condition:`) if its logic is clearly inverted or incorrect. Do not modify it by adding new checker functions or variables.
3.  **Fixing a Flawed Description:** If the accompanying docstring/description is inaccurate or becomes inaccurate after a code fix, you must **REWRITE THE DESCRIPTION** to accurately reflect what the code *actually does*. Rewrite the description to concisely and pricisely describe the function's purpose, based on the code's logic and functionality. The desciption should accurcatly describe the discovered strategy. It should be natural language format " The strategy checks for late stage amide linkage using [AMIDE_COUPLING_REACTION_NAME] in the last two steps of the  synthesis.

You are **STRICTLY FORBIDDEN** from any other form of code editing. Do NOT add new lines, change variable names, or modify the structure of the code in any way other than the allowed improvements.

You will classify each function into one of three categories: **PERFECT**, **GOOD**, or **BAD**, adhering strictly to the following definitions.

---
### Evaluation Criteria & Definitions

**1. PERFECT:**
*   **Code Quality:** The code is flawless, robust, and efficiently uses `checker` functions. Its logic is sound, handles edge cases, and cannot generate false positives.
*   **Strategic Value:** The function identifies a high-level, non-obvious synthetic strategy (e.g., chemoselectivity, late-stage functionalization).
*   **Description:** The description accurately reflects the code's function.

**2. GOOD:**
*   **Code Quality:** The code is robust, correct, and primarily uses `checker` functions. Clever and robust non-checkers for structural features might be okay. Extremely rare edge cases may also be okay.
*   **Strategic Value:** The function identifies a valid but common or lower-level chemical event (e.g., a standard named reaction).
*   **Description:** The description is accurate.

**3. BAD:**
A function is classified as BAD if it meets **ANY** of the following criteria:
*   **Code Quality:** The code is buggy, contains a critical logical flaw, or **fails to use `checker` functions** when they are the superior option.
*   **Strategic Value:** The strategy is trivial, useless, scientifically meaningless, or purely topological.

Consider the following caveats when evaluating:
Low depth corresponds to late-stage in the synthesis, while high depth corresponds to an early-stage in the synthesis. For example, depth 1 is late stage, while depth 6 is early stage.
Reactions are always displayed in the forward direction, so never check for both directions of a reaction. For example, an ester disconnection strategy corresponds to ester formation in the forward direction, not ester cleavage, while ester hydrolysis corresponds to ester cleavage in the forward direction.
For functional group or ring formation checks to be valid, it should either be check_reaction on the whole reaction or check_fg/check_ring on both the reactants and products. Checking only one side is not sufficient.
I am not interested in code which exclusively tests for linear or convergent synthesis, as this is trivial.
---
### Output Format

You MUST respond with a single, valid JSON object and nothing else. The JSON object must have the following structure:

{
  "quality_rating": "string",          // The initial rating: "PERFECT", "GOOD", or "BAD".
  "code_review": "string",             // A brief analysis of the Python code's quality, highlighting any logical flaws.
  "strategy_review": "string",         // An explanation of the chemical strategy's value.
  "suggested_improvement": "string",   // A textual description of WHAT should be improved (code removal, conditional fix, description correction, or a combination).
  "updated_code": "string",            // The full, corrected Python function after applying an allowed code modification. If NO CODE FIX was performed, this MUST be an empty string ("").
  "improved_rating": "string",         // The new quality rating ("PERFECT" or "GOOD") after a code fix. If NO CODE FIX was performed, this MUST be an empty string ("").
  "updated_description": "string"      // A corrected, accurate description. If the original is accurate and remains so, this MUST be an empty string ("").
}
"""
GEMINI_PRO_FILTERING_PROMPT = """You are an expert computational chemist and a senior software engineer with deep knowledge of synthetic organic chemistry. Your task is to analyze, evaluate, and refine Python functions that have been automatically generated to identify specific synthetic strategies from chemical reaction data. The primary aim of these synthetic strategy functions is to annotate a large volume of data and use the resulting one-hot "strategy vectors" for a variety of machine learning tasks. For this end, False Positives and False Negatives will be damaging to the final outcome and MUST be avoided.

Many submitted functions will use a wrapper function (e.g., main(route)) to perform a traversal (e.g., Depth-First Search) over an entire synthesis tree. Your focus for analysis and modification is the inner, recursive function (e.g., dfs_traverse) that processes individual molecules or reactions. The wrapper function itself should not be modified, unless it is to parameterize it with a list of entities created via "Refactoring for Enumeration".

Each step's analysis function will have access to the full reaction object, the current step's depth, and the total number of steps in the synthesis (max_depth). You must ensure this information is correctly propagated and used within the target analysis function.

Evaluation Process & Allowed Improvements
You must meticulously analyze the target function's code, its description, and the chemical strategy it claims to identify. You are permitted to make improvements of the following types, which can occur independently or concurrently.

Propagating Context (Special Initial Step): If the target analysis function (e.g., dfs_traverse(node)) does not have access to reaction, depth, and max_depth, your first modification should be to update its signature and the corresponding call site within the wrapper to pass these parameters correctly. This is a permitted and often necessary refactoring.
Fixing Code by Removal: You can REMOVE any number of lines from the function if they are:
Buggy or Logically Flawed: The code contains clear errors.
Inefficient or Redundant: The code performs unnecessary checks. A common example is checking for a functional group that is already implicitly handled by another, more specific check.
A Source of False Positives: The code uses overly broad or non-specific conditions that are likely to incorrectly flag reactions. A very common source of false positives is chemically incorrect and overly permissive use of FG pattern checks.
Fixing Code by Conditional Modification: You can MODIFY a single conditional statement (if condition:) if its logic is clearly inverted or incorrect.
Constraint: Do not add new checker functions or variables to the condition. The modification should only correct the existing logic using existing elements.
Refactoring for Enumeration (Special Case):
Trigger: This powerful modification is ONLY allowed when the function description specifies a single, specific chemical entity (e.g., "Checks for pyrrole formation") but the code implementation checks for an explicit, well-defined list of related entities (e.g., in a checker call like check_ring(..., ['pyridine', 'pyrrole', 'piperidine'])). This rule does not apply to broad, abstract categories like "any aromatic ring."
Action:
Isolate the List: Move the list of chemical entity strings (e.g., ['pyridine', 'pyrrole', 'piperidine']) outside the function definition, creating a new module-level constant variable (e.g., HETEROCYCLES_OF_INTEREST = [...]).
Update the Code: Modify the function's internal logic to reference this new module-level list instead of the hardcoded one.
Update the Description: Rewrite the description to serve as a template. It should accurately state the general purpose and explicitly reference the list of items being checked. For example: "Checks for the formation of specific heterocyclic rings, including pyridine, pyrrole, and piperidine." The description should exactly match the logic of the code.
Fixing a Flawed Description: If the accompanying docstring/description is inaccurate, or becomes inaccurate after a code fix (including the refactoring above), you must REWRITE THE DESCRIPTION.
Goal: The new description must be concise, precise, and perfectly reflect what the final, corrected code actually does.
You are STRICTLY FORBIDDEN from any form of code editing not explicitly defined in the 'Allowed Improvements' section. Do NOT add new functions, change variable or function names, or alter the fundamental control flow of the code, except as required by the rules above.

Evaluation Criteria & Definitions
You will classify each function into one of three categories: PERFECT, GOOD, or BAD.

---
### Evaluation Criteria & Definitions

**1. PERFECT:**
*   **Code Quality:** The code is flawless, robust, and efficiently uses `checker` functions. Its logic is chemically and computationally sound, handles edge cases, and cannot generate false positives.
*   **Strategic Value:** The function identifies a high-level, non-obvious synthetic strategy (e.g., chemoselectivity, late-stage functionalization).
*   **Description:** The description accurately reflects the code's function.

**2. GOOD:**
*   **Code Quality:** The code is robust, correct, and primarily uses `checker` functions. Clever and robust non-checkers for structural features might be okay. Extremely rare edge cases may also be okay.
*   **Strategic Value:** The function identifies a valid but common or lower-level chemical event (e.g., a standard named reaction).
*   **Description:** The description is accurate.

**3. BAD:**
A function is classified as BAD if it meets **ANY** of the following criteria:
*   **Code Quality:** The code is buggy, contains a critical logical flaw, or **fails to use `checker` functions** when they are the superior option.
*   **Strategic Value:** The strategy is trivial, useless, scientifically meaningless, or purely topological.

Critical Chemical Caveats
Reaction Direction is FORWARD: All reactions are forward synthetic steps.
Synthesis Stages and depth: depth = 1 is the FINAL step (late-stage). depth = max_depth is the FIRST step (early-stage).
Checker Hierarchy: Use of the checker API is strongly preferred over hardcoded SMARTS.
Formation/Cleavage Checks: Must confirm presence/absence on both reactant and product sides.
### Output Format

You MUST respond with a single, valid JSON object and nothing else. The JSON object must have the following structure:

{
  "quality_rating": "string",          // The initial rating: "PERFECT", "GOOD", or "BAD".
  "code_review": "string",             // A brief analysis of the Python code's quality, highlighting any logical flaws.
  "strategy_review": "string",         // An explanation of the chemical strategy's value.
  "suggested_improvement": "string",   // A textual description of WHAT should be improved (code removal, conditional fix, description correction, or a combination).
  "updated_code": "string",            // The full, corrected Python function after applying an allowed code modification. If NO CODE FIX was performed, this MUST be an empty string ("").
  "improved_rating": "string",         // The new quality rating ("PERFECT" or "GOOD") after a code fix. If NO CODE FIX was performed, this MUST be an empty string ("").
  "updated_description": "string"      // A corrected, accurate description. If the original is accurate and remains so, this MUST be an empty string ("").
  "chemist_prompts": "string"          // a few natural language ways a chemist would use to **request** this strategy from a synthesis model - " I want the synthesis to use a ... strategy with {REACTION TYPE} reaction" as an example.
}
"""# (Content omitted for brevity)


def load_function_names(json_file_path: str) -> List[str]:
    """Loads a list of function names from a specific JSON file."""
    try:
        with open(json_file_path, 'r') as f:
            data = json.load(f)
        function_list = list(set(data["selected_functions"]))
        if not isinstance(function_list, list):
             raise TypeError("The value for 'selected_functions' is not a list.")
        print(f"Successfully loaded {len(function_list)} function names from {json_file_path}.")
        return function_list
    except FileNotFoundError:
        print(f"FATAL ERROR: The input file '{json_file_path}' was not found.")
        return []
    except (KeyError, json.JSONDecodeError, TypeError) as e:
        print(f"FATAL ERROR: Could not parse '{json_file_path}'. Details: {e}")
        return []


# MODIFIED: Removed the hardcoded default path. It must now be provided.
def load_function_code(function_name: str, base_path: str) -> Tuple[str, str]:
    """Load the Python file for a given function name from the specified base path."""
    file_path = os.path.join(base_path, f"{function_name}.py")
    if not os.path.exists(file_path):
        # This is now the expected behavior and will be handled in the calling function
        return "", ""
    with open(file_path, 'r') as f:
        code_content = f.read()
    functions = extract_functions(code_content)
    main_function = next((func_code for func_name, func_code in functions if func_name == "main"), None)
    if main_function is None:
        return "", code_content
    return main_function, code_content

# MODIFIED: Accepts model and system_prompt as parameters
def call_openrouter_api(code: str, model_name: str, system_prompt: str) -> str:
    """Call OpenRouter API."""
    try:
        user_prompt = f"Here is the function to analyze:\n\n```python\n{code}\n```"
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={
                "Authorization": f"Bearer {OPENROUTER_API_KEY}",
                "Content-Type": "application/json"
            },
            data=json.dumps({
                "model": model_name,
                "messages": [
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_prompt}
                ],
                "temperature": 0.1,
                "reasoning": {"max_tokens": 8192},
                # This is the key to getting reliable JSON output
                "response_format": {"type": "json_object"}
            }),
            timeout=180
        )
        response.raise_for_status()
        api_response = response.json()
        return api_response['choices'][0]['message']['content']
    except requests.exceptions.RequestException as e:
        print(f"Error calling OpenRouter API: {e}")
        return ""
    except (KeyError, IndexError) as e:
        print(f"Error parsing API response structure: {e}")
        return ""

def parse_api_json_response(response: str) -> Dict[str, Any]:
    """Parses the JSON string response from the API."""
    error_payload = { "quality_rating": "UNKNOWN", "code_review": "Error parsing API response.", "strategy_review": "", "suggested_improvement": "", "updated_code": "", "improved_rating": "", "updated_description": "" }
    if not response:
        return error_payload
    try:
        # The response from chat/completions is often a JSON string, not an object with 'text'
        json_str = response.strip()
        if json_str.startswith("```json"):
            json_match = re.search(r'```json\s*([\s\S]*?)\s*```', json_str)
            if json_match:
                json_str = json_match.group(1)
        parsed_data = json.loads(json_str)
        return parsed_data
    except json.JSONDecodeError:
        print(f"  -> Warning: Could not decode JSON from response. Content: {response[:300]}...")
        return error_payload

# MODIFIED: Accepts all necessary parameters to avoid hardcoding
def process_single_function(function_name: str, source_code_dir: str, model_name: str, system_prompt: str) -> Tuple[str, Dict[str, Any]]:
    """Worker function that processes a single function file."""
    main_code, original_full_code = load_function_code(function_name, base_path=source_code_dir)
    if not main_code:
        print(f"  -> WARNING: Code file not found or 'main' function missing for {function_name} at {os.path.join(source_code_dir, function_name + '.py')}")
        return function_name, {"status": "error", "reason": "Could not load function code."}
    
    response_text = call_openrouter_api(main_code, model_name, system_prompt)
    if not response_text:
        return function_name, {"status": "error", "reason": "API call returned empty response", "original_full_code": original_full_code}
    
    evaluation_data = parse_api_json_response(response_text)
    result = {
        "status": "success",
        "original_full_code": original_full_code,
        "evaluation": evaluation_data
    }
    return function_name, result

# MODIFIED: The main entry point is now fully parameterized.
def process_functions(
    functions_json_path: str,
    source_code_dir: str,
    output_file: str,
    model_name: str,
    system_prompt: str
):
    """Main function to process all functions in parallel and save rated code."""
    function_names = load_function_names(functions_json_path)
    if not function_names:
        print("No function names were loaded. Exiting.")
        return {}
    
    print(f"Starting parallel processing for {len(function_names)} functions with {MAX_CONCURRENT_REQUESTS} workers.")
    print(f"Using model: {model_name}")

    results = {}
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENT_REQUESTS) as executor:
        # Use a lambda to pass the extra, static arguments to the worker function
        future_to_function = {
            executor.submit(
                process_single_function,
                name,
                source_code_dir,
                model_name,
                system_prompt
            ): name for name in function_names
        }
        for future in tqdm.tqdm(as_completed(future_to_function), total=len(function_names), desc="Evaluating Functions"):
            try:
                function_name, result_data = future.result()
                results[function_name] = result_data
            except Exception as exc:
                function_name = future_to_function[future]
                print(f"'{function_name}' generated an exception: {exc}")
                results[function_name] = {"status": "error", "reason": str(exc)}

    print(f"\nAll tasks completed. Saving full evaluation results to {output_file}...")
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    print("Save complete.")
    