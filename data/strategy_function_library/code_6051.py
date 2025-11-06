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


def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthetic route involves formation of aromatic C-N bonds.
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

    c_n_bond_count = 0
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal c_n_bond_count, result, findings_json

        if node.get("type") == "reaction":
            # The original implementation was flawed for manually parsing SMILES.
            # It is replaced by a call to the checker API using the reaction object.
            reaction = node.get("reaction_object")
            if reaction and checker.check_bond_formation(reaction, "[c:1]-[N:2]"):
                c_n_bond_count += 1
                findings_json["atomic_checks"]["named_reactions"].append("aromatic_C-N_bond_formation")

        # Traverse children
        for child in node.get("children", []):
            if node.get("type") == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from root
    dfs_traverse(route)

    if c_n_bond_count >= 1:
        result = True
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "aromatic_C-N_bond_formation",
                "operator": ">=",
                "value": 1
            }
        })

    return result, findings_json
