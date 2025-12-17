from typing import Tuple, Dict, List
import copy
from rdkit.Chem import AllChem, rdFMCS
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


SONOGASHIRA_REACTIONS = [
    "Sonogashira acetylene_aryl halide",
    "Sonogashira alkyne_aryl halide",
    "Sonogashira acetylene_aryl OTf",
    "Sonogashira alkyne_aryl OTf",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects the use of a diaryl alkyne linker, identified by its presence in the final product or by the use of specific Sonogashira coupling reactions to form it. The Sonogashira variants checked are: Sonogashira acetylene_aryl halide, Sonogashira alkyne_aryl halide, Sonogashira acetylene_aryl OTf, and Sonogashira alkyne_aryl OTf.
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

    has_diaryl_alkyne = False
    has_alkyne_forming_reaction = False

    # This is a placeholder for a real calculation of the synthesis depth.
    # The key change is propagating depth information during traversal.
    max_depth = 10 

    def dfs_traverse(node, depth, max_depth):
        nonlocal has_diaryl_alkyne, has_alkyne_forming_reaction, findings_json

        # Check for the linker structure ONLY in the final product (depth=1).
        if node["type"] == "mol" and depth == 1:
            mol_smiles = node["smiles"]
            try:
                # Replaced direct RDKit calls with the preferred checker API.
                if checker.check_substructure("c-C#C-c", mol_smiles):
                    print(f"Detected diaryl alkyne linker in molecule: {mol_smiles}")
                    has_diaryl_alkyne = True
                    findings_json["atomic_checks"]["functional_groups"].append("c-C#C-c")
            except Exception as e:
                print(f"Error processing molecule {mol_smiles}: {e}")

        elif node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rxn_smiles = node["metadata"]["mapped_reaction_smiles"]
            try:
                # Refactored to check against a predefined list for clarity.
                for name in SONOGASHIRA_REACTIONS:
                    if checker.check_reaction(name, rxn_smiles):
                        print(f"Detected Sonogashira coupling reaction: {rxn_smiles}")
                        has_alkyne_forming_reaction = True
                        findings_json["atomic_checks"]["named_reactions"].append(name)
                        break
            except Exception as e:
                print(f"Error processing reaction {rxn_smiles}: {e}")

        # Continue traversing the tree, propagating context.
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            
            dfs_traverse(child, new_depth, max_depth)

    # Start traversal from the root, passing initial context.
    dfs_traverse(route, 1, max_depth)

    strategy_detected = has_diaryl_alkyne or has_alkyne_forming_reaction

    if strategy_detected:
        print("Detected alkyne linker strategy")

    # Populate structural constraints based on the flags
    if has_diaryl_alkyne:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "c-C#C-c",
                "position": "final_product"
            }
        })
    if has_alkyne_forming_reaction:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "sonogashira_coupling",
                "operator": ">=",
                "value": 1
            }
        })

    return strategy_detected, findings_json
