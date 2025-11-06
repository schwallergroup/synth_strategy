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
    Detects the use of silyl protection in a synthetic route.
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

    has_silyl_protection = False

    def dfs_traverse(node, depth=0):
        nonlocal has_silyl_protection, findings_json

        if node["type"] == "reaction":
            rsmi = node.get("metadata", {}).get("rsmi", "")
            if not rsmi:
                return

            # Check for silyl protection reaction
            if checker.check_reaction("Alcohol protection with silyl ethers", rsmi):
                print(f"Detected silyl protection in reaction: {rsmi}")
                has_silyl_protection = True
                findings_json["atomic_checks"]["named_reactions"].append("Alcohol protection with silyl ethers")

            # Alternative check for silyl protective group
            # This catches cases where the specific reaction might not be recognized
            try:
                product_smiles = rsmi.split(">")[-1]
                reactants_smiles = rsmi.split(">")[0].split(".")

                # Check for silyl protective group more generally
                if checker.check_fg("Silyl protective group", product_smiles):
                    # Check if it wasn't present in all reactants (indicating protection occurred)
                    if not all(
                        checker.check_fg("Silyl protective group", r) for r in reactants_smiles
                    ):
                        print(f"Detected silyl protective group formation in reaction: {rsmi}")
                        has_silyl_protection = True
                        findings_json["atomic_checks"]["functional_groups"].append("Silyl protective group")
            except Exception as e:
                print(f"Error checking for silyl groups: {e}")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if has_silyl_protection:
        print("Silyl protection strategy detected")
    else:
        print("No silyl protection detected")

    return has_silyl_protection, findings_json
