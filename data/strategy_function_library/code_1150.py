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
    This function detects if heterocycle formation occurs in the late stage of synthesis.
    It specifically looks for benzoxazole formation in the second half of the synthesis.
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

    benzoxazole_pattern = Chem.MolFromSmarts("c1ccc2c(c1)nco2")
    max_depth = 0
    heterocycle_formation_depth = None
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal max_depth, heterocycle_formation_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Check if benzoxazole is in product but not in reactants
            product_mol = Chem.MolFromSmiles(product_smiles)
            if product_mol and product_mol.HasSubstructMatch(benzoxazole_pattern):
                benzoxazole_in_reactants = False
                for reactant_smiles in reactants_smiles:
                    reactant_mol = Chem.MolFromSmiles(reactant_smiles)
                    if reactant_mol and reactant_mol.HasSubstructMatch(benzoxazole_pattern):
                        benzoxazole_in_reactants = True
                        break

                if not benzoxazole_in_reactants:
                    heterocycle_formation_depth = depth
                    if "benzoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("benzoxazole")
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical)
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if heterocycle formation occurs in the second half of synthesis
    if heterocycle_formation_depth is not None and max_depth > 0:
        relative_position = heterocycle_formation_depth / max_depth
        if relative_position <= 0.5:
            result = True
            findings_json["structural_constraints"].append({
                "type": "positional",
                "details": {
                    "target": "benzoxazole_formation",
                    "position": "late_stage"
                }
            })

    return result, findings_json
