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


PRESERVED_HETEROCYCLES = [
    "benzothiazole",
    "benzofuran",
    "thiazole",
    "benzothiophene",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a specific heterocycle is present in every molecule along the main synthetic path.
    The heterocycles checked are defined in the PRESERVED_HETEROCYCLES list, including:
    benzothiazole, benzofuran, thiazole, and benzothiophene.
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

    heterocycles_by_depth = {}
    result = False

    def dfs_traverse(node, current_depth=0):
        nonlocal findings_json
        if node["type"] == "mol" and "smiles" in node:
            if current_depth not in heterocycles_by_depth:
                heterocycles_by_depth[current_depth] = {
                    hc: False for hc in PRESERVED_HETEROCYCLES
                }
            
            for hc in PRESERVED_HETEROCYCLES:
                if checker.check_ring(hc, node["smiles"]):
                    heterocycles_by_depth[current_depth][hc] = True
                    if hc not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(hc)

        # Determine the depth for the recursive call based on the current node's type
        next_depth = current_depth
        if node["type"] != "reaction": # If current node is 'mol' or other non-reaction type
            next_depth = current_depth + 1

        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    dfs_traverse(route)

    depths = sorted(heterocycles_by_depth.keys())
    if not depths:
        return False, findings_json

    for hc in PRESERVED_HETEROCYCLES:
        # Check if this heterocycle exists at all
        hc_exists = any(
            data.get(hc, False) for data in heterocycles_by_depth.values()
        )

        if hc_exists:
            # Check if it's present at ALL depths where molecules exist
            is_preserved = all(
                heterocycles_by_depth[d].get(hc, False) for d in depths
            )
            if is_preserved:
                result = True
                # Add the structural constraint if preservation is confirmed
                findings_json["structural_constraints"].append(
                    {
                        "type": "preservation",
                        "details": {
                            "target_type": "ring_system",
                            "targets": [
                                "benzothiazole",
                                "benzofuran",
                                "thiazole",
                                "benzothiophene"
                            ],
                            "condition": "any",
                            "scope": "all_stages"
                        }
                    }
                )
                # Since we return True if *any* heterocycle is preserved, we can break here
                break

    return result, findings_json
