import os
import json
import requests
import tqdm
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import List, Optional, Dict, Any, Tuple
import threading
from functools import partial
import copy
import shutil
# For the progress bar (pip install tqdm)

from .llm_utils import extract_functions

OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY")
if not OPENROUTER_API_KEY:
    raise ValueError("OPENROUTER_API_KEY environment variable not set. Please set your OpenRouter API key.")

# --- Concurrency Configuration ---
MAX_CONCURRENT_REQUESTS = 50 # Adjusted for stability

CODE_BASE_PATH = "../data/merged_good_perf_depth_refactored"
# The JSON file containing the structured metadata for each function.
METADATA_DB_PATH = "../data/retrieval/function_metadata_database.json"
# The directory where the new, modified Python files will be saved.
OUTPUT_DIR = "../data/modified_functions_with_structured_output_pro"
# The path to save the summary JSON of the modification run.
RESULTS_JSON_PATH = "modification_run_structured_results.json"

# --- MODIFIED: New prompt for refactoring function return values ---
FUNCTION_MODIFICATION_PROMPT = """
You are an expert Python programmer specializing in refactoring scientific code. Your task is to modify a given Python function to change its return value from a simple boolean to a more informative, structured JSON object representing the findings.

**Context:**
Each function you will modify has a corresponding structured JSON object that defines the "strategy" it is trying to detect. You will be given both the Python function code and its corresponding JSON strategy definition.

**The Goal:**
The input function currently returns a single boolean value (`True`/`False`). You must refactor it to return a tuple containing two elements:
1.  The original boolean value.
2.  A dictionary (`findings_json`) that is a *subset* of the original JSON strategy definition, containing only the elements that were actually detected.

**Instructions:**
1.  **Analyze the inputs:** You will receive the Python code and its JSON strategy definition.
2.  **Modify the function signature:**
    - Add the necessary imports: `from typing import Tuple, Dict, List` and `import copy`.
    - Change the function's return type hint to `-> Tuple[bool, Dict]`.
3.  **Modify the function body:**
    - At the beginning of the function, initialize a `findings_json` dictionary. This should be a deep copy of a predefined empty template that matches the strategy JSON structure.
      ```python
      findings_template = {
          "atomic_checks": {
              "named_reactions": [],
              "ring_systems": [],
              "functional_groups": []
          },
          "structural_constraints": []
      }
      findings_json = copy.deepcopy(findings_template)
      ```
    - Locate every conditional check that results in the main boolean flag being set to `True`.
    - When a check passes, you must now do two things:
        a. Set the main boolean flag to `True`.
        b. Add a record of the successful check to the `findings_json` dictionary.
    - **How to record findings:**
        - **Atomic Checks:** If a specific functional group, reaction, or ring system is found, add its *string name* to the corresponding list within `findings_json["atomic_checks"]`.
            - e.g., `if checker.check_fg("Boronic acid", reactant): ... findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")`
            - e.g., For loops over reactions: `if checker.check_reaction(r_name, rsmi): ... findings_json["atomic_checks"]["named_reactions"].append(r_name)`
        - **Structural Constraints:** If a logical condition representing a structural constraint is met (like `boronic_acid_formed = True`), find the corresponding constraint object in the *original input strategy JSON* and append that entire object to the `findings_json["structural_constraints"]` list.
        - **Positinal Constraint Reformating:** This output must be standardied into the following format    
            ```json
            {
                "type": "positional",
                "targets" : []
                "position": "<position_string>"
            }
            ```
             - <position_string> should be "late_stage", "early_stage", "mid_synthesis", "final_step", "initial_step", "any". Early stage should be used if specified in the function. When handlign operators, remember that low depth is late stage, and high depth is early stage.
        - **Co-Occurrence Constraint Reformating:** This output must be standardied into the following format
            ```json
            {
                "type": "co-occurrence",
                "targets": []
            }
        - **Count Constraint Reformating:** This output must be standardied into the following format
            ```json
            {
                "type": "count",
                "targets": [],
                "operator": ">=",
                "value": 2
            }
          },
        -** Negative Constraint Reformating:** This output must be standardied into the following format
            ```json
            {
                "type": "negation",
                "targets": []
            }
            ```
        - **Other Constraints:** If the original JSON has other types of constraints, you can add them as needed, but ensure they are formatted consistently as above. Prefer categories like "positional", "co-occurrence", "count", and "negation" for clarity.
    These refactoring constraints MUST be followed to allow structured retrieval of function outputs at a later date.
    
4.  **Modify the return statement:**
    - Change the final `return result` to `return result, findings_json`.
5.  **Maintain Code Integrity:** Do NOT change any other part of the function's core logic. The goal is a targeted refactoring of the return value only.

**Example:**
*Original Strategy JSON:*
```json
{
  "atomic_checks": {
    "functional_groups": ["Boronic acid", "Boc"],
    "named_reactions": ["Suzuki coupling", "Grignard reaction", "Boc protection", "Swern oxidation"]
  },
  "constraints": [
    {
      "id": "c1",
      "type": "co-occurrence",
      "details": {
        "targets": ["boronic_acid_formation", "boronic_acid_consumption_in_coupling"]
      }
    },
    {
      "id": "c2",
      "type": "positional",
      "details": {
        "target": "Boc_protection",
        "position": {
          "operator": "<=",
          "value": 2 
        }
      }
    },
    {
        "id": "c3",
        "type": "count",
        "details": {
            "target": "Grignard_reaction",
            "operator": ">=",
            "value": 2
        }
    },
    {
        "id": "c4",
        "type": "negation",
        "details": {
            "target": "Swern_oxidation"
        }
    }
  ]
}
```
Example Refactored Function:


def main(route: Dict[str, Any], strategy_json: Dict[str, Any]) -> Tuple[bool, Dict]:
    # --- Initialization ---
    findings_template = {
        "atomic_checks": {"named_reactions": [], "ring_systems": [], "functional_groups": []},
        "constraints": [] # Changed from "structural_constraints" to be more general
    }
    findings_json = copy.deepcopy(findings_template)

    # --- Logic flags to track conditions ---
    boronic_acid_formed = False
    boronic_acid_used_in_coupling = False
    boc_protection_in_early_stage = False
    grignard_reaction_count = 0
    swern_oxidation_found = False

    # --- Traversal Logic (Conceptual) ---
    # This part remains conceptually the same, it just sets more flags.
    def dfs_traverse(node: Dict[str, Any], depth: int):
        nonlocal boronic_acid_formed, boronic_acid_used_in_coupling, findings_json
        nonlocal boc_protection_in_early_stage, grignard_reaction_count, swern_oxidation_found

        if node.get("type") == "reaction":
            rsmi = node.get("rsmi", "")
            
            # Atomic checks (unchanged logic)
            if checker.check_fg("Boronic acid", node):
                boronic_acid_formed = True
                findings_json["atomic_checks"]["functional_groups"].append("Boronic acid")
            
            if checker.check_reaction("Suzuki coupling", rsmi):
                 boronic_acid_used_in_coupling = True
                 findings_json["atomic_checks"]["named_reactions"].append("Suzuki coupling")

            # Logic for new constraints
            if checker.check_reaction("Boc protection", rsmi):
                findings_json["atomic_checks"]["named_reactions"].append("Boc protection")
                # Check if the position constraint (<= 2 steps from the start) is met
                if depth <= 2: # This simulates the 'early_stage' condition
                    boc_protection_in_early_stage = True

            if checker.check_reaction("Grignard reaction", rsmi):
                 grignard_reaction_count += 1
                 findings_json["atomic_checks"]["named_reactions"].append("Grignard reaction")

            if checker.check_reaction("Swern oxidation", rsmi):
                swern_oxidation_found = True
                findings_json["atomic_checks"]["named_reactions"].append("Swern oxidation")


        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route, 0)

    # --- Post-Traversal Constraint Evaluation and Refactoring ---
    # This section iterates through the *original* constraints, checks if the
    # logic was satisfied, and appends the *reformatted* finding.

    all_constraints_met = True
    
    for constraint in strategy_json.get("constraints", []):
        constraint_met = False
        
        # 1. Structural Constraint (Co-occurrence)
        if constraint["type"] == "co-occurrence":
            if boronic_acid_formed and boronic_acid_used_in_coupling:
                constraint_met = True
                # Per the new rule, find the original object and append it.
                # NOTE: The prompt is slightly ambiguous. For maximum utility, we reformat
                # ALL constraint findings into the new standard format.
                findings_json["constraints"].append({
                    "type": "co-occurrence",
                    "targets": constraint["details"]["targets"]
                })

        # 2. Positional Constraint
        elif constraint["type"] == "positional":
            if boc_protection_in_early_stage:
                constraint_met = True
                # Reformat according to the specified positional structure
                position_str = "early_stage" # Mapped from operator '<='
                findings_json["constraints"].append({
                    "type": "positional",
                    "targets": [constraint["details"]["target"]],
                    "position": position_str
                })
        
        # 3. Count Constraint
        elif constraint["type"] == "count":
            op = constraint["details"]["operator"]
            val = constraint["details"]["value"]
            if (op == ">=" and grignard_reaction_count >= val):
                constraint_met = True
                # Reformat according to the specified count structure
                findings_json["constraints"].append({
                    "type": "count",
                    "targets": [constraint["details"]["target"]],
                    "operator": op,
                    "value": val
                })

        # 4. Negative Constraint
        elif constraint["type"] == "negation":
            if not swern_oxidation_found:
                constraint_met = True
                # Reformat according to the specified negation structure
                findings_json["constraints"].append({
                    "type": "negation",
                    "targets": [constraint["details"]["target"]]
                })
        
        if not constraint_met:
            all_constraints_met = False

    # --- Final Return Statement ---
    # The final boolean result depends on all constraints being satisfied.
    return all_constraints_met, findings_json

###
Expected Findings JSON Output:
```json
        {
    "atomic_checks": {
        "named_reactions": ["Suzuki coupling", "Boc protection", "Grignard reaction"],
        "ring_systems": [],
        "functional_groups": ["Boronic acid"]
    },
    "constraints": [
        {
        "type": "co-occurrence",
        "targets": ["boronic_acid_formation", "boronic_acid_consumption_in_coupling"]
        },
        {
        "type": "positional",
        "targets": ["Boc_protection"],
        "position": "early_stage"
        },
        {
        "type": "count",
        "targets": ["Grignard_reaction"],
        "operator": ">=",
        "value": 2
        },
        {
        "type": "negation",
        "targets": ["Swern_oxidation"]
        }
    ]
    }
```
Output Format:
You MUST respond with ONLY a single JSON object. Do not include any preamble, explanations, or markdown formatting. The JSON object must have the following structure:
{
"modification_successful": boolean,
"reason": "A brief explanation of why the modification was or was not successful.",
"modified_code": "The complete, refactored Python code for the new function, as a single string. This should be an empty string if modification_successful is false."
}
"""

# --- Helper Functions ---

def load_metadata_database(db_path: str) -> Dict[str, Any]:
        """Loads the entire function metadata JSON file into a dictionary."""
        print(f"Loading function metadata from: {db_path}")
        if not os.path.exists(db_path):
            print(f"FATAL ERROR: Metadata database not found at {db_path}")
            return {}
        try:
            with open(db_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
                # Assuming the data is a list of objects, convert it to a dict keyed by function_id
                metadata_db = {item['function_id']: item for item in data}
                print(f"Successfully loaded metadata for {len(metadata_db)} functions.")
                return metadata_db
        except (json.JSONDecodeError, KeyError) as e:
            print(f"FATAL ERROR: Could not parse metadata database {db_path}: {e}")
            return {}

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

def call_llm_api(code: str, metadata_json: str) -> str:
    """Calls the specified LLM API with the modification prompt, code, and metadata."""
    if not OPENROUTER_API_KEY or OPENROUTER_API_KEY == "YOUR_API_KEY_HERE":
        print("FATAL ERROR: OPENROUTER_API_KEY is not set.")
        return ""
    try:
        user_prompt = f"{FUNCTION_MODIFICATION_PROMPT} Here is the function to analyze:\n\npython\n{code}\n. \nAnd here is the metadata JSON:\n\n{metadata_json}\n\n"
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

def parse_api_json_response(response: str) -> Dict[str, Any]:
    """Parses the JSON string response from the API."""
    error_payload = {"modification_successful": False, "reason": "Error parsing API response.", "modified_code": ""}
    if not response: return error_payload
    try:
        return json.loads(response)
    except json.JSONDecodeError:
        print(f"  -> Warning: Could not decode JSON from response. Content: {response[:300]}...")
        return error_payload

# --- Main Processing Logic ---

def process_single_function(function_name: str, code_base_path: str, metadata_db: Dict[str, Any]) -> Tuple[str, Dict[str, Any]]:
    """
    Processes a single function: loads code and metadata, sends to LLM, and parses the result.
    """
    print(f"Processing '{function_name}'...")

    # Step 1: Load code
    full_code = load_function_code(function_name, code_base_path)
    if not full_code:
        return function_name, {"status": "error", "reason": "Could not load code file."}

    # Step 2: Load corresponding metadata
    metadata = metadata_db.get(function_name)
    if not metadata:
        return function_name, {"status": "error", "reason": f"No metadata found for function_id '{function_name}'."}
    metadata_str = json.dumps(metadata, indent=2)

    # Step 3: Call LLM with both code and metadata
    response_text = call_llm_api(full_code, metadata_str)
    if not response_text:
        return function_name, {"status": "error", "reason": "API call failed or returned empty."}

    # Step 4: Parse response
    modification_data = parse_api_json_response(response_text)
    return function_name, {
        "status": "success",
        "original_function_name": function_name,
        "modification_details": modification_data
    }

def process_and_modify_functions(code_dir: str, metadata_path: str, output_dir: str, output_json_path: str):
    """
    Main function to process all functions in parallel, save modified code, and copy non-modified code.
    """
    metadata_db = load_metadata_database(metadata_path)
    if not metadata_db:
        print("Exiting due to metadata loading failure.")
        return

    function_names = list(metadata_db.keys())
    print(f"Will process {len(function_names)} functions based on metadata.")

    # Use the passed `output_dir`
    os.makedirs(output_dir, exist_ok=True)
    print(f"Output function files will be saved in '{output_dir}/'")

    results = {}
    worker_func = partial(process_single_function, code_base_path=code_dir, metadata_db=metadata_db)
    
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENT_REQUESTS) as executor:
        future_to_function = {executor.submit(worker_func, name): name for name in function_names}
        for future in tqdm.tqdm(as_completed(future_to_function), total=len(function_names), desc="Modifying Functions with LLM"):
            try:
                name, result_data = future.result()
                results[name] = result_data
            except Exception as exc:
                name = future_to_function[future]
                print(f"'{name}' generated an exception: {exc}")
                results[name] = {"status": "error", "reason": str(exc)}

    print(f"\nAll tasks completed. Saving summary to {output_json_path}...")
    with open(output_json_path, 'w') as f:
        json.dump(results, f, indent=2)
    print("Summary saved.")

    print("\nSaving output function files...")
    modified_count, not_modified_count, error_count, total_files_written = 0, 0, 0, 0

    for original_name, data in results.items():
        if data.get("status") == "error":
            error_count += 1
            print(f"  -> SKIPPING '{original_name}' due to processing error: {data.get('reason')}")
            continue

        details = data.get("modification_details", {})
        
        if details.get("modification_successful"):
            modified_count += 1
            new_code = details.get("modified_code")
            if new_code:
                # Use the passed `output_dir`
                save_path = os.path.join(output_dir, f"{original_name}.py")
                with open(save_path, 'w', encoding='utf-8') as f:
                    f.write(new_code)
                total_files_written += 1
            else:
                error_count += 1
        else:
            not_modified_count += 1
            original_code_path = os.path.join(code_dir, f"{original_name}.py")
            if os.path.exists(original_code_path):
                # Use the passed `output_dir`
                output_path = os.path.join(output_dir, f"{original_name}.py")
                try:
                    shutil.copy(original_code_path, output_path)
                    total_files_written += 1
                except Exception as e:
                    error_count += 1
            else:
                error_count += 1

    print("\n--- Final Summary ---")
    print(f"Functions processed:        {len(results)}")
    print(f"Successfully modified:      {modified_count}")
    print(f"Copied (not modified):    {not_modified_count}")
    print(f"Errors during processing:   {error_count}")
    print(f"Total files written:        {total_files_written} in '{output_dir}'")