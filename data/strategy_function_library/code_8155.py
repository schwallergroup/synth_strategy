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
    This function detects if the synthesis uses an alkyne as a rigid linker
    between aromatic systems in the final product.
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

    has_alkyne_linker = False

    def dfs_traverse(node, depth=0):
        nonlocal has_alkyne_linker, findings_json

        if node["type"] == "mol" and depth == 0:  # Final product
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check for alkyne connecting two aromatic systems
                alkyne_linker_pattern = Chem.MolFromSmarts("[c,n][C]#[C][c,n]")
                if mol.HasSubstructMatch(alkyne_linker_pattern):
                    has_alkyne_linker = True
                    findings_json["atomic_checks"]["functional_groups"].append("alkyne")
                    findings_json["structural_constraints"].append({
                        "type": "positional",
                        "details": {
                            "target": "alkyne_linker_between_aromatics",
                            "position": "last_stage"
                        }
                    })
                    print("Alkyne linker detected in final product")

        # Continue traversing
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # chemical or mol
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return has_alkyne_linker, findings_json