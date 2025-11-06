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
    Detects a strategy where a nitro group is preserved throughout the synthesis
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

    # Track reactions where nitro group is preserved
    preserved_nitro_reactions = []

    def dfs_traverse(node, depth=0):
        nonlocal findings_json
        if node["type"] == "reaction":
            try:
                # Extract reactants and product from reaction SMILES
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1].split(".")

                # Check for nitro group in products and reactants
                product_has_nitro = any(
                    checker.check_fg("Nitro group", smi) for smi in product_smiles if smi
                )
                reactants_has_nitro = any(
                    checker.check_fg("Nitro group", smi) for smi in reactants_smiles if smi
                )

                # If nitro is preserved in this reaction, record it
                if product_has_nitro and reactants_has_nitro:
                    preserved_nitro_reactions.append(depth)
                    if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
            except Exception:
                pass

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'chemical' or other types
                dfs_traverse(child, depth + 1)

    # Check if target molecule has nitro group
    target_has_nitro = False
    if route["type"] == "mol":
        target_has_nitro = checker.check_fg("Nitro group", route["smiles"])
        if target_has_nitro:
            if "Nitro group" not in findings_json["atomic_checks"]["functional_groups"]:
                findings_json["atomic_checks"]["functional_groups"].append("Nitro group")
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "Nitro group in final product"
                    ]
                }
            })
        print(f"Target molecule has nitro group: {target_has_nitro}")

    # Start traversal
    dfs_traverse(route)

    # Strategy is present if:
    # 1. Target molecule has nitro group
    # 2. At least one reaction preserves the nitro group
    # 3. The preservation spans multiple reaction steps (at least 2 depths)
    result = False
    if target_has_nitro and preserved_nitro_reactions:
        # Check if preservation spans multiple depths
        unique_depths = set(preserved_nitro_reactions)
        if len(unique_depths) >= 2:
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "unique_depths_of_nitro_group_preservation",
                    "operator": ">=",
                    "value": 2
                }
            })
            print(
                f"Strategy detected: Nitro group preservation throughout synthesis at depths {sorted(unique_depths)}"
            )
            result = True

    return result, findings_json