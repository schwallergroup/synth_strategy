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
    This function detects if the synthesis includes sequential deprotection of Boc and Methoxy protecting groups.
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

    deprotection_depths = []
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal deprotection_depths, findings_json

        if node["type"] == "reaction":
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for Boc deprotection
            boc_pattern = Chem.MolFromSmarts("[#6]([#6])([#6])([#6])[#8][#6](=[#8])[#7]")
            has_boc_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(boc_pattern):
                    has_boc_in_reactants = True
                    findings_json["atomic_checks"]["functional_groups"].append("Boc group")
                    break

            product_mol = Chem.MolFromSmiles(product)
            has_boc_in_product = product_mol and product_mol.HasSubstructMatch(boc_pattern)

            if has_boc_in_reactants and not has_boc_in_product:
                deprotection_depths.append((depth, "Boc"))
                findings_json["atomic_checks"]["named_reactions"].append("Boc deprotection")

            # Check for methoxy deprotection
            methoxy_pattern = Chem.MolFromSmarts("[#6][#8][#6]:[#6]")

            has_methoxy_in_reactants = False
            for reactant in reactants:
                reactant_mol = Chem.MolFromSmiles(reactant)
                if reactant_mol and reactant_mol.HasSubstructMatch(methoxy_pattern):
                    has_methoxy_in_reactants = True
                    findings_json["atomic_checks"]["functional_groups"].append("aryl methyl ether")
                    break

            has_methoxy_in_product = product_mol and product_mol.HasSubstructMatch(methoxy_pattern)

            if has_methoxy_in_reactants and not has_methoxy_in_product:
                deprotection_depths.append((depth, "Methoxy"))
                findings_json["atomic_checks"]["named_reactions"].append("Methoxy deprotection")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction":  # If current node is chemical, increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    # Sort deprotection events by depth
    deprotection_depths.sort()

    # Check if we have at least 2 different deprotection events
    if len(deprotection_depths) >= 2:
        # Check if they are different types
        deprotection_types = set([d[1] for d in deprotection_depths])
        if len(deprotection_types) >= 2:
            print(f"Detected sequential deprotection strategy: {deprotection_depths}")
            result = True
            # Add structural constraint if both Boc and Methoxy deprotections are found
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "Boc deprotection",
                        "Methoxy deprotection"
                    ]
                }
            })

    return result, findings_json