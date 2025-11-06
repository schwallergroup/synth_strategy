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

HETEROCYCLES_OF_INTEREST = [
    "indole",
    "quinoline",
    "isoquinoline",
    "pyridine",
    "pyrimidine",
    "pyrazole",
    "imidazole",
    "thiazole",
    "oxazole",
    "triazole",
    "tetrazole",
    "furan",
    "thiophene",
    "pyrrole",
    "morpholine",
    "piperidine",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the final product contains multiple heterocyclic systems,
    specifically looking for indole and quinolone/quinoline cores.
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

    has_multiple_heterocycles = False

    def dfs_traverse(node, depth=0):
        nonlocal has_multiple_heterocycles, findings_json

        if node["type"] == "mol" and depth == 0:  # Final product
            try:
                total_heterocycles = 0
                for ring_type in HETEROCYCLES_OF_INTEREST:
                    ring_indices = checker.get_ring_atom_indices(ring_type, node["smiles"])
                    if ring_indices:
                        # Add to findings_json if found
                        if ring_type not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append(ring_type)
                    total_heterocycles += len(ring_indices)

                if total_heterocycles >= 2:
                    has_multiple_heterocycles = True
                    # Add structural constraint if condition met
                    findings_json["structural_constraints"].append({
                        "type": "count",
                        "details": {
                            "target": "heterocycle",
                            "operator": ">=",
                            "value": 2
                        }
                    })
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "heterocycle_count_check",
                            "position": "last_stage"
                        }
                    })

            except Exception as e:
                print(f"Error analyzing heterocycles: {e}")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    dfs_traverse(route)
    return has_multiple_heterocycles, findings_json
