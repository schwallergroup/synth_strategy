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


HETEROCYCLE_TYPES = [
    "furan",
    "pyran",
    "pyrrole",
    "pyridine",
    "pyrazole",
    "imidazole",
    "oxazole",
    "thiazole",
    "pyrimidine",
    "pyrazine",
    "pyridazine",
    "triazole",
    "tetrazole",
    "pyrrolidine",
    "piperidine",
    "piperazine",
    "morpholine",
    "thiomorpholine",
    "indole",
    "quinoline",
    "isoquinoline",
    "thiophene",
    "benzothiophene",
    "benzoxazole",
    "benzothiazole",
    "benzimidazole",
]

def main(route) -> Tuple[bool, Dict]:
    """
    Detects a synthesis that involves at least 3 different types of heterocycles from a predefined list. This function traverses the entire synthetic route, inspects each molecule for the presence of heterocycles specified in the `HETEROCYCLE_TYPES` list, and returns True if the diversity threshold is met.
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

    # Track heterocycles found
    heterocycles_found = set()

    def dfs_traverse(node, depth=0):
        nonlocal heterocycles_found, findings_json

        if node["type"] == "mol" and "smiles" in node:
            mol_smiles = node["smiles"]

            # Check for each heterocycle type
            for heterocycle in HETEROCYCLE_TYPES:
                if checker.check_ring(heterocycle, mol_smiles):
                    heterocycles_found.add(heterocycle)
                    if heterocycle not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append(heterocycle)
                    print(f"Found heterocycle: {heterocycle} in molecule {mol_smiles}")

        # Traverse children
        for child in node.get("children", []):
            # Determine the new depth based on the current node's type
            new_depth = depth if node["type"] == "reaction" else depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    result = False
    # Return True if we found at least 3 different heterocycles
    if len(heterocycles_found) >= 3:
        print(f"Found heterocycle-rich synthesis with: {', '.join(heterocycles_found)}")
        result = True
        findings_json["structural_constraints"].append({"type": "count", "details": {"target": "heterocycle_diversity", "operator": ">=", "value": 3}})
    else:
        print(
            f"Found only {len(heterocycles_found)} heterocycles: {', '.join(heterocycles_found) if heterocycles_found else 'none'}"
        )
    
    return result, findings_json