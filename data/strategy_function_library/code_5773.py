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
from steerable_retro.utils.check import Check
from steerable_retro.utils import fuzzy_dict, check

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
    This function detects if a benzothiophene ring system is formed
    in any step of the synthesis prior to the final step.
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

    # Track if we found heterocycle formation
    found_formation = False

    def dfs_traverse(node, depth=0):
        nonlocal found_formation, findings_json

        # Non-final stage reactions (depth >= 2)
        if node["type"] == "reaction" and depth >= 2:
            # Extract reactants and product
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if product contains benzothiophene
            if checker.check_ring("benzothiophene", product_smiles):
                findings_json["atomic_checks"]["ring_systems"].append("benzothiophene")

                # Check if any reactant already has the complete benzothiophene
                reactant_has_pattern = False
                for reactant in reactants_smiles:
                    if checker.check_ring("benzothiophene", reactant):
                        reactant_has_pattern = True
                        break

                # If no reactant has benzothiophene but product does, we found a formation
                if not reactant_has_pattern:
                    found_formation = True
                    # Add the structural constraint if the formation is detected in a non-final stage
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "benzothiophene_ring_formation",
                            "position": "not_last_stage"
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            # New logic: depth increases only when going from chemical to reaction
            # Depth remains same when going from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    print("Starting traversal to find non-final step benzothiophene formation")
    dfs_traverse(route)
    print(f"Non-final step benzothiophene formation detected: {found_formation}")
    return found_formation, findings_json
