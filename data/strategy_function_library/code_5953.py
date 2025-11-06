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
    Detects preservation of nitro group throughout the synthesis
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

    nitro_group_present_in_product_count = 0
    total_reactions = 0

    def dfs_traverse(node, depth=0):
        nonlocal total_reactions, nitro_group_present_in_product_count, findings_json

        if node["type"] == "reaction":
            total_reactions += 1
            if "mapped_reaction_smiles" in node["metadata"]:
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                product = rsmi.split(">")[-1]

                # Check if product has nitro group
                product_mol = Chem.MolFromSmiles(product)
                nitro_pattern = Chem.MolFromSmarts("[N+](=O)[O-]")

                if product_mol and product_mol.HasSubstructMatch(nitro_pattern):
                    nitro_group_present_in_product_count += 1
                    if "nitro" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("nitro")
                    print(f"Detected nitro group at depth {depth}")

        # Continue traversing
        for child in node.get("children", []):
            # New logic: depth increases only when traversing from chemical to reaction
            # Depth remains the same when traversing from reaction to chemical
            new_depth = depth + 1 if node["type"] != "reaction" else depth
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if nitro group is present in all reactions
    result = nitro_group_present_in_product_count == total_reactions and total_reactions > 0

    if total_reactions > 0:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "any_reaction",
                "operator": ">",
                "value": 0
            }
        })

    if nitro_group_present_in_product_count == total_reactions:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_without_nitro_in_product",
                "operator": "==",
                "value": 0
            }
        })

    print(f"Nitro group present in {nitro_group_present_in_product_count}/{total_reactions} reactions")
    print(f"Nitro group preservation strategy detected: {result}")
    return result, findings_json