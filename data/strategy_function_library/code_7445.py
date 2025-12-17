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

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route uses an early Boc protection
    strategy (depth >= 2).
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

    boc_protection_early = False

    def dfs_traverse(node, depth=0):
        nonlocal boc_protection_early, findings_json

        print(f"Traversing node at depth {depth}, type: {node['type']}")

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"].get("rsmi", "")
                if not rsmi:
                    print(f"No reaction SMILES found at depth {depth}")
                    return

                print(f"Checking reaction at depth {depth}: {rsmi[:50]}...")

                # Extract reactants and product
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check if any reactant already has Boc
                reactants_with_boc = False
                for reactant in reactants:
                    if reactant.strip() and checker.check_fg("Boc", reactant):
                        reactants_with_boc = True
                        findings_json["atomic_checks"]["functional_groups"].append("Boc")
                        break

                # Check if product has Boc
                product_has_boc = checker.check_fg("Boc", product)
                if product_has_boc:
                    findings_json["atomic_checks"]["functional_groups"].append("Boc")

                print(
                    f"Reactants with Boc: {reactants_with_boc}, Product has Boc: {product_has_boc}"
                )

                # If Boc was added (not in reactants but in product) and it's at an early stage
                if product_has_boc and not reactants_with_boc and depth >= 2:
                    print(f"Detected early Boc protection at depth {depth}")
                    boc_protection_early = True
                    # Record the structural constraint
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "Boc protection",
                            "position": "early_stage (depth >= 2)"
                        }
                    })
                    return
            except Exception as e:
                print(f"Error in Boc protection detection: {e}")

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction":  # If current node is 'chemical'
            next_depth = depth + 1

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)
    print(f"Final result: {boc_protection_early}")
    return boc_protection_early, findings_json
