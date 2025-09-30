# ### File: scripts/run_pipeline.py (CORRECTED) ###

import os
import sys
import json
import argparse
import shutil
from typing import List

# --- Add project root to sys.path ---
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, PROJECT_ROOT)

try:
    from src.synth_strategy.llm import function_filtering
    from src.synth_strategy.llm import function_database_extraction
    from src.synth_strategy.llm import function_return_modification
except ImportError as e:
    print(f"Error: Could not import pipeline modules. Ensure correct project structure. Error: {e}")
    sys.exit(1)

# ... (parse_arguments, generate_initial_function_list, create_passed_functions_list remain the same) ...
def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments for the pipeline."""
    parser = argparse.ArgumentParser(description="A 4-stage pipeline for processing synthetic strategy functions.")
    
    parser.add_argument(
        "--initial-functions-json",
        type=str,
        default="data/filtering_data/all_functions.json",
        help="Path to the initial JSON file containing function names. If it doesn't exist, it will be generated automatically from the source code directory."
    )
    parser.add_argument(
        "--source-code-dir",
        type=str,
        required=True,
        help="Path to the directory containing the original Python files for all functions."
    )
    parser.add_argument(
        "--base-output-dir",
        type=str,
        default="./pipeline_output",
        help="Base directory where processed code artifacts for each stage will be stored."
    )
    parser.add_argument(
        "--reports-dir",
        type=str,
        default="data/filtering_data",
        help="Directory to store all intermediate JSON reports (e.g., filter ratings, modification summaries)."
    )
    parser.add_argument(
        "--run-stages",
        type=int,
        nargs='+',
        default=[1, 2, 3, 4],
        help="A list of integers specifying which stages to run (e.g., --run-stages 1 3 4)."
    )
    parser.add_argument(
        "--flash-filter-level",
        type=str,
        choices=["perfect", "good_and_perfect"],
        default="good_and_perfect",
        help="For Stage 1 (Gemini Flash), specifies which function ratings to keep."
    )
    return parser.parse_args()

def generate_initial_function_list(source_dir: str, output_json_path: str):
    """Scans a directory for .py files and creates the initial function list JSON."""
    print(f"Initial function list not found at '{output_json_path}'. Generating it now...")
    if not os.path.isdir(source_dir):
        print(f"FATAL ERROR: The specified source code directory does not exist: {source_dir}")
        sys.exit(1)
        
    function_names = [os.path.splitext(f)[0] for f in os.listdir(source_dir) if f.endswith('.py')]
    
    if not function_names:
        print(f"FATAL ERROR: No Python files found in {source_dir}. Cannot generate function list.")
        sys.exit(1)

    os.makedirs(os.path.dirname(output_json_path), exist_ok=True)
    
    with open(output_json_path, 'w') as f:
        json.dump({"selected_functions": function_names}, f, indent=2)
        
    print(f"Successfully generated list of {len(function_names)} functions at '{output_json_path}'.")

def create_passed_functions_list(
    report_path: str,
    source_code_dir: str,
    output_dir: str,
    ratings_to_keep: List[str]
) -> tuple[str, str]:
    """
    Reads a filtering report, copies passing files from the source_code_dir to a new
    directory, and creates a new JSON list for the next stage.
    """
    if not os.path.exists(report_path):
        print(f"  ERROR: Report file not found at {report_path}. Cannot proceed.")
        return "", ""

    with open(report_path, 'r') as f:
        results = json.load(f)

    passed_function_names = []
    for func_name, data in results.items():
        if data.get("status") != "success": continue
        eval_data = data.get("evaluation", {})
        final_rating = eval_data.get("improved_rating") or eval_data.get("quality_rating")
        
        if final_rating and final_rating.upper() in ratings_to_keep:
            passed_function_names.append(func_name)

    print(f"  Found {len(passed_function_names)} functions with ratings {ratings_to_keep} to pass to the next stage.")

    new_code_dir = os.path.join(output_dir, "functions")
    os.makedirs(new_code_dir, exist_ok=True)
    for func_name in passed_function_names:
        src_file = os.path.join(source_code_dir, f"{func_name}.py")
        dst_file = os.path.join(new_code_dir, f"{func_name}.py")
        if os.path.exists(src_file):
            shutil.copy(src_file, dst_file)
        else:
            print(f"  WARNING: Source code for '{func_name}' not found at {src_file}")

    new_functions_json_path = os.path.join(output_dir, "passed_functions.json")
    with open(new_functions_json_path, 'w') as f:
        json.dump({"selected_functions": passed_function_names}, f, indent=2)
        
    return new_functions_json_path, new_code_dir


def main():
    args = parse_arguments()
    
    os.makedirs(args.base_output_dir, exist_ok=True)
    os.makedirs(args.reports_dir, exist_ok=True)
    print(f"Pipeline code artifacts will be saved in: {os.path.abspath(args.base_output_dir)}")
    print(f"Pipeline JSON reports will be saved in:   {os.path.abspath(args.reports_dir)}")

    if not os.path.exists(args.initial_functions_json):
        generate_initial_function_list(args.source_code_dir, args.initial_functions_json)
        
    current_functions_json = args.initial_functions_json
    current_source_code_dir = args.source_code_dir

    # =========================================================================
    # STAGE 1: Filtering with Gemini Flash
    # =========================================================================
    if 1 in args.run_stages:
        print("\n" + "="*50)
        print("--- RUNNING STAGE 1: Filtering with Gemini Flash ---")
        print("="*50)
        
        stage1_output_dir = os.path.join(args.base_output_dir, "1_flash_filter_output")
        os.makedirs(stage1_output_dir, exist_ok=True)
        stage1_report_path = os.path.join(args.reports_dir, "flash_filter_report.json")
        
        print(f"Input functions: {current_functions_json}")
        print(f"Input source code from: {current_source_code_dir}")
        print(f"Output report to: {stage1_report_path}")

        # MODIFIED: Call the parameterized function directly. No more monkey-patching.
        function_filtering.process_functions(
            functions_json_path=current_functions_json,
            source_code_dir=current_source_code_dir,
            output_file=stage1_report_path,
            model_name="google/gemini-2.5-flash", # Or your preferred flash model
            system_prompt=function_filtering.GEMINI_FLASH_FILTERING_PROMPT
        )
        
        ratings_to_keep = ["PERFECT"]
        if args.flash_filter_level == "good_and_perfect":
            ratings_to_keep.append("GOOD")

        current_functions_json, current_source_code_dir = create_passed_functions_list(
            report_path=stage1_report_path,
            source_code_dir=args.source_code_dir,
            output_dir=stage1_output_dir,
            ratings_to_keep=ratings_to_keep
        )
        print("--- Stage 1 Complete ---")

    # =========================================================================
    # STAGE 2: Filtering with Gemini Pro
    # =========================================================================
    if 2 in args.run_stages:
        print("\n" + "="*50)
        print("--- RUNNING STAGE 2: Filtering with Gemini Pro ---")
        print("="*50)

        if not current_functions_json or not os.path.exists(current_functions_json):
            print("ERROR: No functions passed Stage 1. Skipping Stage 2.")
        else:
            stage2_output_dir = os.path.join(args.base_output_dir, "2_pro_filter_output")
            os.makedirs(stage2_output_dir, exist_ok=True)
            stage2_report_path = os.path.join(args.reports_dir, "pro_filter_report.json")
            
            print(f"Input functions: {current_functions_json}")
            print(f"Input source code from: {current_source_code_dir}")
            print(f"Output report to: {stage2_report_path}")
            
            # MODIFIED: Call the parameterized function directly.
            function_filtering.process_functions(
                functions_json_path=current_functions_json,
                source_code_dir=current_source_code_dir,
                output_file=stage2_report_path,
                model_name="google/gemini-2.5-pro", # Or your preferred pro model
                system_prompt=function_filtering.GEMINI_PRO_FILTERING_PROMPT
            )
            
            current_functions_json, current_source_code_dir = create_passed_functions_list(
                report_path=stage2_report_path,
                source_code_dir=current_source_code_dir,
                output_dir=stage2_output_dir,
                ratings_to_keep=["PERFECT", "GOOD"]
            )
        print("--- Stage 2 Complete ---")

    # ... (Stages 3 and 4 remain unchanged as they were already correctly parameterized) ...
    # =========================================================================
    # STAGE 3: Metadata Generation
    # =========================================================================
    if 3 in args.run_stages:
        print("\n" + "="*50)
        print("--- RUNNING STAGE 3: Metadata Generation ---")
        print("="*50)

        if not current_source_code_dir or not os.path.exists(current_source_code_dir) or not os.listdir(current_source_code_dir):
            print("ERROR: No function files available from previous stages. Skipping Stage 3.")
        else:
            stage3_output_dir = os.path.join(args.base_output_dir, "3_metadata_output")
            os.makedirs(stage3_output_dir, exist_ok=True)
            stage3_db_path = os.path.join(stage3_output_dir, "function_metadata_database.json")
            
            print(f"Input code directory: {current_source_code_dir}")
            print(f"Output database: {stage3_db_path}")

            function_database_extraction.generate_function_metadata_db(
                code_dir=current_source_code_dir,
                output_json_path=stage3_db_path
            )
        print("--- Stage 3 Complete ---")

    # =========================================================================
    # STAGE 4: Return Value Modification
    # =========================================================================
    if 4 in args.run_stages:
        print("\n" + "="*50)
        print("--- RUNNING STAGE 4: Return Value Modification ---")
        print("="*50)
        
        stage3_db_path = os.path.join(args.base_output_dir, "3_metadata_output", "function_metadata_database.json")

        if not os.path.exists(stage3_db_path):
             print(f"ERROR: Metadata database not found at {stage3_db_path}. Run Stage 3 first. Skipping Stage 4.")
        elif not current_source_code_dir or not os.path.exists(current_source_code_dir) or not os.listdir(current_source_code_dir):
            print("ERROR: No function files available from previous stages. Skipping Stage 4.")
        else:
            stage4_output_dir = os.path.join(args.base_output_dir, "4_modified_functions_output")
            os.makedirs(stage4_output_dir, exist_ok=True)
            stage4_report_path = os.path.join(args.reports_dir, "modification_report.json")
            
            function_return_modification.OUTPUT_DIR = stage4_output_dir
            
            print(f"Input code directory: {current_source_code_dir}")
            print(f"Input metadata DB: {stage3_db_path}")
            print(f"Output directory for modified code: {stage4_output_dir}")
            print(f"Output report to: {stage4_report_path}")
            
            function_return_modification.process_and_modify_functions(
                code_dir=current_source_code_dir,
                output_dir=stage4_output_dir,
                metadata_path=stage3_db_path,
                output_json_path=stage4_report_path
            )
        print("--- Stage 4 Complete ---")

    print("\nPIPELINE FINISHED.")

if __name__ == "__main__":
    main()