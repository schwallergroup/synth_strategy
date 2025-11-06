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


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if any starting material in the synthesis contains a trifluoromethyl group.
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

    has_trifluoromethyl = False

    def dfs_traverse(node, depth=0):
        nonlocal has_trifluoromethyl, findings_json

        if node["type"] == "mol" and "smiles" in node and not node.get("children"):
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C]([F])([F])[F]")):
                has_trifluoromethyl = True
                findings_json["atomic_checks"]["functional_groups"].append("trifluoromethyl")
                findings_json["structural_constraints"].append({"type": "positional", "details": {"target": "trifluoromethyl", "position": "starting_material"}})

        # Traverse children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    return has_trifluoromethyl, findings_json