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
    This function detects benzyl protection/deprotection as a synthetic strategy.
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

    benzyl_deprotection_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal benzyl_deprotection_detected, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for benzyl group in reactant but not in product (deprotection)
                reactant_mol = Chem.MolFromSmiles(reactants[0])
                product_mol = Chem.MolFromSmiles(product)

                if reactant_mol and product_mol:
                    benzyl_pattern = Chem.MolFromSmarts("[#8][C][c]1[cH][cH][cH][cH][cH]1")

                    if (
                        reactant_mol.HasSubstructMatch(benzyl_pattern)
                        and not product_mol.HasSubstructMatch(benzyl_pattern)
                    ):
                        benzyl_deprotection_detected = True
                        findings_json["atomic_checks"]["functional_groups"].append("benzyl_ether_or_ester")
                        findings_json["atomic_checks"]["named_reactions"].append("benzyl_deprotection")
                        # Assuming benzene ring is part of benzyl, if detected, add it
                        benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
                        if reactant_mol.HasSubstructMatch(benzene_pattern):
                            findings_json["atomic_checks"]["ring_systems"].append("benzene")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from a reaction node
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from a chemical node
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if benzyl_deprotection_detected:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "benzyl_deprotection",
                "operator": ">=",
                "value": 1
            }
        })

    return benzyl_deprotection_detected, findings_json
