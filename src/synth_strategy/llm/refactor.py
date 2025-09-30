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

from .llm_utils import extract_functions

# --- OpenRouter API Configuration ---
OPENROUTER_API_KEY = os.getenv("OPENROUTER_API_KEY")
if not OPENROUTER_API_KEY:
    raise ValueError("OPENROUTER_API_KEY environment variable not set. Please set your OpenRouter API key.")

# --- Concurrency and Rate Limit Configuration ---
MAX_CONCURRENT_REQUESTS = 50
MAX_RPD = 50000
# --- Directory and File Configuration ---
# MODIFIED: Input and Output directories updated as requested
INPUT_DIR = "../data/merged_good_perf"
OUTPUT_DIR = "../data/merged_good_perf_depth_refactored"
RESULTS_JSON_PATH = "merged_good_perf_depth_refactored.json"


# --- LLM Prompt Configuration ---
# MODIFIED: Prompt is now targeted for the specific depth refactoring task
REFACTORING_PROMPT = """
Your task is to rewrite the provided Python script to change its depth calculation logic.

**Instructions:**
1.  Locate the `dfs_traverse` function within the script.
2.  Modify its recursive call. The new rule is that depth should **only increase** when traversing from a chemical node to a reaction node. The depth should **remain the same** when traversing from a reaction node to a chemical node.
    - **Incorrect (Old Logic):** `dfs_traverse(child, depth + 1)` for all children.
    - **Correct (New Logic):** The depth passed to the recursive call should be `depth` if the current node's type is 'reaction', and `depth + 1` if the current node's type is not 'reaction' (e.g., 'chemical').
3.  Do not change any other part of the code's logic. Ensure all original imports and helper code outside the `dfs_traverse` function are preserved.

Before you start, please read the entire script to understand its structure and logic. Your goal is to ensure that the depth calculation reflects the new rule while keeping the rest of the code intact.

**Output Format:**
You MUST respond with ONLY a single JSON object. Do not include any preamble, explanations, or markdown formatting around the JSON block. The JSON object must have this exact structure:

{
  "refactored_code": "The complete, rewritten Python code as a single string."
}
"""

# --- Helper Functions ---

def load_function_names_from_directory(directory_path: str) -> List[str]:
    """Scans a directory for .py files and returns their base names."""
    print(f"Scanning directory for function names: {directory_path}")
    if not os.path.isdir(directory_path):
        print(f"FATAL ERROR: The specified directory does not exist: {directory_path}")
        return []

    function_names = [
        os.path.splitext(f)[0]
        for f in os.listdir(directory_path)
        if f.endswith('.py') and os.path.isfile(os.path.join(directory_path, f))
    ]
    print(f"Found {len(function_names)} files to process.")
    return function_names

def load_function_code(function_name: str, base_path: str) -> Optional[str]:
    """Loads the full Python file content for a given function name."""
    file_path = os.path.join(base_path, f"{function_name}.py")
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            return f.read()
    except Exception as e:
        print(f"  -> ERROR: Failed to read file {file_path}: {e}")
        return None

def call_llm_api(code: str) -> str:
    """
    Calls the specified LLM API with the refactoring prompt.
    NOTE: The API call itself is unchanged, as requested.
    """
    if not OPENROUTER_API_KEY or OPENROUTER_API_KEY == "YOUR_API_KEY_HERE":
        print("FATAL ERROR: OPENROUTER_API_KEY is not set.")
        return ""
    try:
        user_prompt = f"{REFACTORING_PROMPT}\n\nHere is the script to refactor:\n```python\n{code}\n```"
        response = requests.post(
            url="https://openrouter.ai/api/v1/chat/completions",
            headers={"Authorization": f"Bearer {OPENROUTER_API_KEY}"},
            data=json.dumps({
                "model": "google/gemini-2.5-flash",
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
    """Parses the JSON string response from the API for the refactored code."""
    error_payload = {"refactored_code": None, "reason": "Error parsing API response."}
    if not response: return error_payload
    try:
        data = json.loads(response)
        if "refactored_code" not in data:
            data["reason"] = "JSON response missing 'refactored_code' key."
        return data
    except json.JSONDecodeError:
        print(f"  -> Warning: Could not decode JSON from response. Content: {response[:300]}...")
        return error_payload

# --- Main Processing Logic ---

def process_single_file(function_name: str, code_base_path: str) -> Tuple[str, Dict[str, Any]]:
    """Processes a single file: loads its code, sends it to the LLM, and parses the result."""
    # print(f"Processing '{function_name}'...")
    full_code = load_function_code(function_name, code_base_path)
    if not full_code:
        return function_name, {"status": "error", "reason": "Could not load code file."}

    response_text = call_llm_api(full_code)
    if not response_text:
        return function_name, {"status": "error", "reason": "API call failed or returned empty."}

    parsed_data = parse_api_json_response(response_text)
    return function_name, {
        "status": "success" if parsed_data.get("refactored_code") else "error",
        "original_function_name": function_name,
        "details": parsed_data
    }

def process_and_refactor_files(input_dir: str, output_dir: str, output_json_path: str):
    """
    Main function to process all files in parallel, save refactored code, and create a summary.
    """
    function_names = load_function_names_from_directory(input_dir)
    if not function_names:
        print("No files to process. Exiting.")
        return

    os.makedirs(output_dir, exist_ok=True)
    print(f"Output files will be saved in '{output_dir}/'")

    results = {}
    with ThreadPoolExecutor(max_workers=MAX_CONCURRENT_REQUESTS) as executor:
        future_to_function = {
            executor.submit(process_single_file, name, input_dir): name for name in function_names
        }
        for future in tqdm.tqdm(as_completed(future_to_function), total=len(function_names), desc="Refactoring Files with LLM"):
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

    print("\nSaving refactored output files...")
    refactored_count, error_count = 0, 0

    for original_name, data in results.items():
        if data.get("status") == "success":
            details = data.get("details", {})
            new_code = details.get("refactored_code")

            if new_code:
                save_path = os.path.join(output_dir, f"{original_name}.py")
                try:
                    with open(save_path, 'w', encoding='utf-8') as f:
                        f.write(new_code)
                    refactored_count += 1
                except IOError as e:
                    print(f"  -> ERROR: Could not write file {save_path}: {e}")
                    error_count += 1
            else:
                 print(f"  -> ERROR: '{original_name}' had success status but no code.")
                 error_count += 1
        else:
            print(f"  -> SKIPPING '{original_name}' due to error: {data.get('reason') or data.get('details', {}).get('reason')}")
            error_count += 1

    print("\n--- Final Summary ---")
    print(f"Total files processed:      {len(results)}")
    print(f"Successfully refactored:    {refactored_count}")
    print(f"Errors during processing:   {error_count}")
    print(f"-----------------------------")
    print(f"Total files written to '{output_dir}': {refactored_count}")
    print("-----------------------------")

if __name__ == "__main__":
    process_and_refactor_files(
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        output_json_path=RESULTS_JSON_PATH
    )
    print("\nScript finished.")