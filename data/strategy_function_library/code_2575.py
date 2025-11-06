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
    Detects if a beta-lactam moiety is present in the final product and also appears in at least one precursor molecule in the synthesis tree.
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

    # Track beta-lactam presence
    has_beta_lactam_in_final = False
    has_beta_lactam_in_intermediates = False

    # SMARTS pattern for beta-lactam
    beta_lactam_pattern = "[#6]1[#6][#7][#6](=[#8])1"

    def dfs_traverse(node, depth=0):
        nonlocal has_beta_lactam_in_final, has_beta_lactam_in_intermediates, findings_json

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol and mol.HasSubstructMatch(Chem.MolFromSmarts(beta_lactam_pattern)):
                findings_json["atomic_checks"]["functional_groups"].append("beta-lactam")
                if depth == 0:  # Final product
                    has_beta_lactam_in_final = True
                else:  # Intermediate
                    has_beta_lactam_in_intermediates = True

        # Determine the next depth based on the current node's type
        next_depth = depth
        if node["type"] != "reaction": # Depth increases only when going from chemical to reaction
            next_depth = depth + 1

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, next_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if the strategy is present
    preserves_beta_lactam = has_beta_lactam_in_final and has_beta_lactam_in_intermediates

    if has_beta_lactam_in_final:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "beta-lactam",
                "position": "last_stage"
            }
        })
    if has_beta_lactam_in_intermediates:
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "beta-lactam",
                "position": "not_last_stage"
            }
        })

    if preserves_beta_lactam:
        print("Detected beta-lactam preservation strategy")

    return preserves_beta_lactam, findings_json