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
    Detects if no new stereocenters are created during the synthesis. This is confirmed by checking that the number of stereocenters in any intermediate molecule is less than or equal to the number of stereocenters in the final product. A temporary loss of stereocenters is considered acceptable.
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
    result = True

    # Get the target molecule's stereocenter count
    target_mol = Chem.MolFromSmiles(route["smiles"])
    Chem.AssignStereochemistry(target_mol)
    target_stereocenters = len(Chem.FindMolChiralCenters(target_mol, includeUnassigned=False))

    # If target has no stereocenters, there's nothing to preserve
    if target_stereocenters == 0:
        result = False
        return result, findings_json
    else:
        # Add the structural constraint if target has stereocenters
        findings_json["structural_constraints"].append(
            {
                "type": "count",
                "details": {
                    "target": "stereocenters_in_final_product",
                    "operator": ">",
                    "value": 0
                }
            }
        )

    # Track intermediate molecules (non-starting materials)
    intermediates = []

    def dfs_traverse(node, depth=0):
        nonlocal intermediates
        if node["type"] == "mol":
            # Skip starting materials
            if not node.get("in_stock", False):
                mol = Chem.MolFromSmiles(node["smiles"])
                if mol:
                    # Count stereocenters
                    Chem.AssignStereochemistry(mol)
                    chiral_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=False))
                    intermediates.append(chiral_centers)

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else: # Assuming 'mol' or other types
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Check if any intermediate has more stereocenters than the target
    # This would indicate new stereocenters were created
    for count in intermediates:
        if count > target_stereocenters:
            result = False
            return result, findings_json

    # If the loop completes, no intermediate has more stereocenters than the target.
    # This means stereochemistry was preserved or reduced, which is acceptable.
    # Add the structural constraint if the condition holds true for all intermediates
    findings_json["structural_constraints"].append(
        {
            "type": "count",
            "details": {
                "target": "stereocenters_in_intermediate",
                "operator": "<=",
                "value": "stereocenters_in_final_product"
            }
        }
    )
    result = True
    return result, findings_json
