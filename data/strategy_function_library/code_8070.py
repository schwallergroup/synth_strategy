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
    Detects routes where a nitrile group is present in the final product and is largely maintained throughout the synthesis. A route is flagged if: 1. The final target molecule contains a nitrile group. 2. Of all reaction steps that produce a product containing a nitrile, the nitrile was also present in at least one reactant for 75% or more of those steps.
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

    # Track reactions where nitrile is maintained
    reactions_with_nitrile_maintained = 0
    total_reactions_with_nitrile_product = 0

    # Check if the final target molecule contains a nitrile
    target_has_nitrile = checker.check_fg("Nitrile", route["smiles"])
    if target_has_nitrile:
        findings_json["atomic_checks"]["functional_groups"].append("Nitrile")

    def dfs_traverse(node, depth=0):
        nonlocal reactions_with_nitrile_maintained, total_reactions_with_nitrile_product, findings_json

        if node["type"] == "reaction":
            try:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Check if product has nitrile
                product_has_nitrile = checker.check_fg("Nitrile", product_smiles)

                if product_has_nitrile:
                    total_reactions_with_nitrile_product += 1

                    # Check if any reactant has nitrile
                    reactants_have_nitrile = any(
                        checker.check_fg("Nitrile", reactant) for reactant in reactants_smiles
                    )

                    # Nitrile is maintained if it's in both product and at least one reactant
                    if reactants_have_nitrile:
                        reactions_with_nitrile_maintained += 1
            except Exception:
                # Silently ignore errors in malformed reaction data
                pass

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for its children (which are chemicals)
                dfs_traverse(child, depth)
            else:
                # If current node is a chemical, depth increases for its children (which are reactions)
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Return True if:
    # 1. The target molecule has a nitrile group
    # 2. At least 75% of reactions that produce a nitrile maintain it from a reactant
    nitrile_maintenance_ratio = (
        reactions_with_nitrile_maintained / total_reactions_with_nitrile_product
        if total_reactions_with_nitrile_product > 0
        else 1.0
    )

    result = target_has_nitrile and nitrile_maintenance_ratio >= 0.75

    if target_has_nitrile:
        findings_json["structural_constraints"].append(
            {
                "type": "count",
                "details": {
                    "target": "Nitrile_functional_group_on_final_product",
                    "operator": ">=",
                    "value": 1
                }
            }
        )
    
    if nitrile_maintenance_ratio >= 0.75:
        findings_json["structural_constraints"].append(
            {
                "type": "count",
                "details": {
                    "target": "nitrile_maintenance_ratio",
                    "operator": ">=",
                    "value": 0.75
                }
            }
        )

    return result, findings_json
