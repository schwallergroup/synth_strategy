import json
import os
from pathlib import Path
from collections import defaultdict

def extract_strategies_to_files(json_data, output_dir, root_data="/home/dparm/synth_strategy/data"):
    """
    Extract passing code blocks from JSON and write to individual Python files.
    
    Args:
        json_data: The loaded JSON data containing strategies and code blocks
        output_dir: Directory where Python files will be written
        root_data: Root directory for data patterns (default: /home/dparm/synth_strategy/data)
    """
    # Create output directory if it doesn't exist
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Track function names to handle duplicates
    function_name_counts = defaultdict(int)
    
    # Header template for each file
    header = f'''import re
from rdkit import Chem
from collections import deque
import fuzzy_dict
import check

# Setup data paths and resources
root_data = "{root_data}"
fg_args = {{
    "file_path": f"{{root_data}}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}}
reaction_class_args = {{
    "file_path": f"{{root_data}}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}}
ring_smiles_args = {{
    "file_path": f"{{root_data}}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)
checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)

'''
    
    total_strategies = 0
    total_written = 0
    total_skipped = 0
    
    # Iterate through all dictionaries in the JSON data
    for item in json_data:
        # Each item is a dict with numeric string keys like "0", "1", etc.
        for idx_key, strategy_data in item.items():
            print(f"\n=== Processing strategy set '{idx_key}' ===")
            
            code_blocks = strategy_data.get("code_blocks", [])
            final_passed = strategy_data.get("final_passed", [])
            
            # Process each code block
            for code_idx, (code, passed) in enumerate(zip(code_blocks, final_passed)):
                total_strategies += 1
                
                if passed:
                    # Extract function name from code
                    func_name = None
                    for line in code.split('\n'):
                        if line.strip().startswith('def strategy_'):
                            func_name = line.split('(')[0].replace('def ', '').strip()
                            break
                    
                    if not func_name:
                        print(f"  ⚠ Warning: Could not extract function name from code block {code_idx}")
                        func_name = f"strategy_{idx_key}_{code_idx}"
                    
                    # Handle duplicate function names
                    base_func_name = func_name
                    function_name_counts[base_func_name] += 1
                    
                    if function_name_counts[base_func_name] > 1:
                        # Append a suffix to make it unique
                        suffix = function_name_counts[base_func_name] - 1
                        func_name = f"{base_func_name}_v{suffix}"
                        print(f"  ℹ Duplicate function name detected. Renaming to: {func_name}")
                    
                    # Remove duplicate imports from code block
                    code_lines = code.split('\n')
                    filtered_lines = []
                    for line in code_lines:
                        # Skip import lines that are already in header
                        if line.strip().startswith('import re'):
                            continue
                        elif line.strip().startswith('from rdkit import Chem'):
                            continue
                        elif line.strip().startswith('from rdkit.Chem import AllChem'):
                            continue
                        elif line.strip().startswith('from collections import deque'):
                            continue
                        else:
                            filtered_lines.append(line)
                    
                    cleaned_code = '\n'.join(filtered_lines)
                    
                    # Create complete file content
                    file_content = header + cleaned_code
                    
                    # Write to file
                    filename = f"{func_name}.py"
                    filepath = output_path / filename
                    
                    with open(filepath, 'w') as f:
                        f.write(file_content)
                    
                    print(f"  ✓ Written: {filename}")
                    total_written += 1
                else:
                    print(f"  ✗ Skipped code block {code_idx} in set '{idx_key}' (did not pass)")
                    total_skipped += 1
    
    print(f"\n{'='*60}")
    print(f"Extraction Summary:")
    print(f"  Total strategies found: {total_strategies}")
    print(f"  Successfully written: {total_written}")
    print(f"  Skipped (failed): {total_skipped}")
    print(f"{'='*60}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract strategy functions from JSON to individual Python files')
    parser.add_argument('json_file', type=str, help='Path to the JSON file containing strategies')
    parser.add_argument('--output-dir', '-o', type=str, default='data/test',
                        help='Output directory for strategy files (default: ./strategy_function_library)')
    parser.add_argument('--root-data', '-r', type=str, default='/home/dparm/synth_strategy/data',
                        help='Root data directory path (default: /home/dparm/synth_strategy/data)')
    
    args = parser.parse_args()
    
    # Load the JSON data
    with open(args.json_file, 'r') as f:
        data = json.load(f)
    
    # Extract strategies to files
    extract_strategies_to_files(data, args.output_dir, args.root_data)
    
    print(f"\nExtraction complete! Files written to: {args.output_dir}")