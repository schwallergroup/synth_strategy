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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects a synthetic strategy where the final step
    involves reduction of a nitro group to an amine.
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

    final_step_has_nitro_reduction = False

    def dfs_traverse(node, depth=0):
        nonlocal final_step_has_nitro_reduction, findings_json

        # Check if this is the final product (a molecule with no children)
        if node["type"] == "mol" and depth == 0:
            print(f"Found final product molecule: {node['smiles']}")
            # The final product should have a parent reaction
            if node.get("children") and len(node["children"]) > 0:
                child = node["children"][0]  # Get the reaction that produced this molecule
                if child["type"] == "reaction" and "rsmi" in child.get("metadata", {}):
                    rsmi = child["metadata"]["rsmi"]

                    print(f"Analyzing final reaction: {rsmi}")

                    # Check if this is a nitro reduction reaction
                    if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                        print("Confirmed nitro reduction reaction type")
                        final_step_has_nitro_reduction = True
                        findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                        findings_json["structural_constraints"].append({
                            "type": "positional",
                            "details": {
                                "target": "Reduction of nitro groups to amines",
                                "position": "last_stage"
                            }
                        })
                        return

        # If this is a reaction node at depth 0, it might be the final step
        elif node["type"] == "reaction" and depth == 0:
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]

                print(f"Analyzing potential final reaction: {rsmi}")

                # Check if this is a nitro reduction reaction
                if checker.check_reaction("Reduction of nitro groups to amines", rsmi):
                    print("Confirmed nitro reduction reaction type")
                    final_step_has_nitro_reduction = True
                    findings_json["atomic_checks"]["named_reactions"].append("Reduction of nitro groups to amines")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Reduction of nitro groups to amines",
                            "position": "last_stage"
                        }
                    })
                    return

        # Continue traversal
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not a reaction (i.e., it's a chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from root
    dfs_traverse(route)
    print(f"Final result: {final_step_has_nitro_reduction}")

    return final_step_has_nitro_reduction, findings_json
