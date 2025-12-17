from typing import Tuple, Dict, List
import copy
from collections import deque
import rdkit.Chem as Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdChemReactions
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import rdkit.Chem.rdFMCS
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import AllChem, Descriptors, Lipinski
from rdkit.Chem import rdmolops
import re
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem, Descriptors
import traceback
import rdkit
from collections import Counter
from synth_strategy.utils.check import Check
from synth_strategy.utils import fuzzy_dict, check

from pathlib import Path
root_data = Path(__file__).parent.parent

fg_args = {
    "file_path": f"{root_data}/patterns/functional_groups.json",
    "value_field": "pattern",
    "key_field": "name",
}
reaction_class_args = {
    "file_path": f"{root_data}/patterns/smirks.json",
    "value_field": "smirks",
    "key_field": "name",
}
ring_smiles_args = {
    "file_path": f"{root_data}/patterns/chemical_rings_smiles.json",
    "value_field": "smiles",
    "key_field": "name",
}
functional_groups = fuzzy_dict.FuzzyDict.from_json(**fg_args)
reaction_classes = fuzzy_dict.FuzzyDict.from_json(**reaction_class_args)
ring_smiles = fuzzy_dict.FuzzyDict.from_json(**ring_smiles_args)

checker = check.Check(
    fg_dict=functional_groups, reaction_dict=reaction_classes, ring_dict=ring_smiles
)


AMIDATION_REACTIONS = [
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_N",
    "Acylation of Nitrogen Nucleophiles by Acyl/Thioacyl/Carbamoyl Halides and Analogs_OS",
    "Carboxylic acid with primary amine to amide",
    "Ester with primary amine to amide",
    "Acyl chloride with primary amine to amide (Schotten-Baumann)",
    "Acyl chloride with secondary amine to amide",
    "Ester with secondary amine to amide",
    "Ester with ammonia to amide",
    "Acyl chloride with ammonia to amide",
    "Schotten-Baumann_amide",
    "Acylation of primary amines",
    "Acylation of secondary amines",
    "Carboxylic acid to amide conversion",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthetic strategy where an amide is formed in the early stages of synthesis, defined as the earliest two-thirds of the total steps. Amide formation is identified by checking if the reaction matches any in a predefined list of common amidation reaction types, such as 'Carboxylic acid with primary amine to amide' and 'Acyl chloride with primary amine to amide (Schotten-Baumann)'.
    """
    findings_template = {
        "atomic_checks": {
            "named_reactions": [],
            "ring_systems": [],
            "functional_groups": []
        },
        "structural_constraints": []
    }
    findings_json = copy.deepcopy(findings_template)

    amidation_found = False
    amidation_depth = float("inf")
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amidation_found, amidation_depth, max_depth, findings_json

        # Track maximum depth to determine total synthesis length
        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                print(f"Checking reaction at depth {depth}: {rsmi}")

                # Check for amidation reactions using the checker functions
                is_amidation = False
                for reaction_name in AMIDATION_REACTIONS:
                    if checker.check_reaction(reaction_name, rsmi):
                        is_amidation = True
                        if reaction_name not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append(reaction_name)
                        break

                if is_amidation:
                    amidation_found = True
                    amidation_depth = min(amidation_depth, depth)
                    print(f"Amidation found at depth {depth}")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Continue traversing the tree
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            if node["type"] == "reaction":
                # From reaction to chemical, depth remains the same
                dfs_traverse(child, depth)
            else:
                # From chemical to reaction, depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    print(
        f"Max depth: {max_depth}, Amidation depth: {amidation_depth if amidation_found else 'not found'}"
    )

    result = False
    # Consider it early-stage if amidation occurs in the first two-thirds of the synthesis
    if amidation_found and max_depth > 0 and amidation_depth >= (max_depth / 3):
        print(f"Early-stage amidation confirmed: depth {amidation_depth} out of {max_depth}")
        result = True
        # Add structural constraint finding
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "any_amidation_reaction",
                "position": "early_stage",
                "description": "The first occurrence of any amidation reaction must be found in the first two-thirds of the synthesis. This is checked by the condition `amidation_depth >= max_depth / 3`, where `amidation_depth` is the depth of the first amidation and `max_depth` is the total length of the route (depth 0 is the last step)."
            }
        })

    if amidation_found:
        print(f"Amidation found but not in early stage: depth {amidation_depth} out of {max_depth}")
    else:
        print("No amidation reactions found in the route")

    return result, findings_json
