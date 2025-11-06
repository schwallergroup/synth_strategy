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


MAINTAINED_HETEROCYCLES = ['indazole', 'benzotriazole', 'pyridine']

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if a synthesis route maintains at least one of a specific list of
    heterocyclic scaffolds, defined in `MAINTAINED_HETEROCYCLES`, throughout the synthesis.
    The list includes indazole, benzotriazole, and pyridine.
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

    heterocycle_data = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal heterocycle_data, findings_json

        if node["type"] == "mol" and "smiles" in node:
            smiles = node["smiles"]
            try:
                mol_data = {"depth": depth, "smiles": smiles}
                for hetero in MAINTAINED_HETEROCYCLES:
                    if checker.check_ring(hetero, smiles):
                        mol_data[hetero] = True
                        if hetero not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(hetero)
                    else:
                        mol_data[hetero] = False
                heterocycle_data.append(mol_data)
            except Exception:
                pass

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node.get("type") != "reaction": # This covers 'mol' and any other non-reaction type
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Filter only molecule data (this is now redundant but harmless)
    molecule_data = [d for d in heterocycle_data if "type" not in d or d["type"] != "reaction"]

    if len(molecule_data) < 2:
        return False, findings_json

    # Sort by depth
    molecule_data.sort(key=lambda x: x["depth"])

    # Check if heterocycles are maintained throughout
    threshold = max(0.7, 1.0 - 0.05 * len(molecule_data))

    # Check if any heterocycle is consistently present and in the final product
    final_product = molecule_data[0]
    for hetero in MAINTAINED_HETEROCYCLES:
        num_present = sum(1 for m in molecule_data if m.get(hetero, False))
        
        heterocycle_in_final_product = final_product.get(hetero, False)
        
        if len(molecule_data) > 0:
            rate = num_present / len(molecule_data)
        else:
            rate = 0

        heterocycle_high_presence_rate = (rate >= threshold)

        if heterocycle_high_presence_rate and heterocycle_in_final_product:
            result = True
            # Add structural constraints to findings_json
            # Positional constraint
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": [
                        "indazole",
                        "benzotriazole",
                        "pyridine"
                    ],
                    "position": "last_stage",
                    "note": "A heterocycle from the target list must be present in the final product (depth=0)."
                }
            })
            # Count constraint
            findings_json["structural_constraints"].append({
                "type": "count",
                "details": {
                    "target": "a_single_heterocycle_from_list",
                    "operator": ">=",
                    "value": "dynamic_threshold",
                    "scope": "all_molecules",
                    "note": "The rate of presence for a single heterocycle across all molecules must exceed a dynamic threshold, which is at least 0.7 (rate = num_present / num_molecules)."
                }
            })
            # Co-occurrence constraint
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "heterocycle_in_final_product",
                        "heterocycle_high_presence_rate"
                    ],
                    "note": "The positional (presence in final product) and count (high presence rate) constraints must be met for the same heterocycle from the target list."
                }
            })
            return result, findings_json

    return result, findings_json