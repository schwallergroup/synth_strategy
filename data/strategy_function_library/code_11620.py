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
    This function detects a synthetic strategy involving the formation and subsequent
    opening of an epoxide as key intermediate steps.
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

    epoxide_formed = False
    epoxide_opened = False
    epoxide_formation_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal epoxide_formed, epoxide_opened, epoxide_formation_depth, findings_json

        if node["type"] == "reaction":
            if "mapped_reaction_smiles" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0]
                products_smiles = rsmi.split(">")[-1]

                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles.split(".") if r]
                products = [Chem.MolFromSmiles(p) for p in products_smiles.split(".") if p]

                if not all(reactants) or not all(products):
                    return

                # Check for epoxide formation
                epoxide_pattern = Chem.MolFromSmarts("[C]1O[C]1")
                if not any(mol.HasSubstructMatch(epoxide_pattern) for mol in reactants) and any(
                    mol.HasSubstructMatch(epoxide_pattern) for mol in products
                ):
                    epoxide_formed = True
                    epoxide_formation_depth = depth
                    if "epoxide" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("epoxide")
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

                # Check for epoxide opening
                if any(mol.HasSubstructMatch(epoxide_pattern) for mol in reactants) and not any(
                    mol.HasSubstructMatch(epoxide_pattern) for mol in products
                ):
                    epoxide_opened = True
                    if "epoxide" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("epoxide")
                    if "ring_destruction" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_destruction")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is not a reaction (e.g., chemical), increase depth
                new_depth = depth + 1
            dfs_traverse(child, new_depth)

    dfs_traverse(route)

    result = epoxide_formed and epoxide_opened
    if result:
        findings_json["structural_constraints"].append({
            "type": "co-occurrence",
            "details": {
                "targets": [
                    "epoxide_formation",
                    "epoxide_opening"
                ]
            }
        })

    return result, findings_json
