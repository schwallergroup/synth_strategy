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


import checker

def main(route) -> Tuple[bool, Dict]:
    """
    This function detects if the synthesis involves formation of a urea linker
    connecting two aromatic systems.
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

    urea_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal urea_formation_detected, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            DIARYL_UREA_SMARTS = "[c,n][NH][C](=[O])[NH][c,n]"

            # The manual RDKit logic is replaced by a single, robust checker call.
            # This assumes the checker API can operate directly on a reaction SMILES string.
            if checker.check_substructure_formation_from_smiles(rsmi, DIARYL_UREA_SMARTS):
                urea_formation_detected = True
                findings_json["atomic_checks"]["functional_groups"].append("urea")
                findings_json["atomic_checks"]["named_reactions"].append("urea_formation")
                print(f"Urea linker formation detected at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node to a chemical node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node to a reaction node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return urea_formation_detected, findings_json