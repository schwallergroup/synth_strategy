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
    Detects if the synthetic route incorporates a fluorinated group,
    specifically a difluoromethoxy group.
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

    fluorinated_group_found = False

    def dfs_traverse(node, depth=0):
        nonlocal fluorinated_group_found, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"].get("rsmi", "")
            if not rsmi:
                return

            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                # Check for difluoromethoxy group
                difluoromethoxy_pattern = Chem.MolFromSmarts("[O][C]([F])[F]")
                if product and product.HasSubstructMatch(difluoromethoxy_pattern):
                    # Verify it's a formation by checking reactants
                    reactants_have_pattern = any(
                        r and r.HasSubstructMatch(difluoromethoxy_pattern) for r in reactants if r
                    )
                    if not reactants_have_pattern:
                        fluorinated_group_found = True
                        findings_json["atomic_checks"]["functional_groups"].append("difluoromethoxy")
                        # Conceptually, this is a 'difluoromethoxy_formation' reaction
                        findings_json["atomic_checks"]["named_reactions"].append("difluoromethoxy_formation")

            except Exception as e:
                # Error handling should be done by the calling framework
                pass

        # Continue traversing children
        for child in node.get("children", []):
            new_depth = depth
            if node["type"] != "reaction": # Only increase depth if current node is not 'reaction'
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Add structural constraint if the difluoromethoxy group was found
    if fluorinated_group_found:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "difluoromethoxy_formation",
                "operator": ">=",
                "value": 1
            }
        })

    return fluorinated_group_found, findings_json