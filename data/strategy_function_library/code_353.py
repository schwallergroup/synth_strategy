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
    This function detects if a gem-difluoroalkyl group (a carbon bonded to two fluorine atoms) is present in the final product and also present in at least one of its precursors (intermediate or starting material).
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

    difluoro_present_in_final = False
    difluoro_present_in_intermediates = False
    result = False

    def dfs_traverse(node, depth):
        nonlocal difluoro_present_in_final, difluoro_present_in_intermediates, findings_json

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                difluoro_pattern = Chem.MolFromSmarts("[#6]([F])([F])")
                if mol.HasSubstructMatch(difluoro_pattern):
                    findings_json["atomic_checks"]["functional_groups"].append("gem-difluoroalkyl group")
                    if depth == 1:
                        # The root node (depth=1) is the final product
                        difluoro_present_in_final = True
                    else:
                        # All other nodes are intermediates or starting materials
                        difluoro_present_in_intermediates = True

        # Continue traversing, incrementing depth based on node type
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node
                dfs_traverse(child, depth + 1)

    # Start traversal from the root (final product) at depth 1
    dfs_traverse(route, 1)

    # The strategy is present if difluoro is in both final product and at least one precursor
    if difluoro_present_in_final and difluoro_present_in_intermediates:
        result = True
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "gem-difluoroalkyl group in final product",
                    "gem-difluoroalkyl group in precursor"
                ]
            }
        })

    return result, findings_json