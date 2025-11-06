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
    This function detects if a synthetic route employs a strategy where a nitrile group
    is maintained through multiple steps and then converted to an isoxazole in the final step.
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

    nitrile_pattern = Chem.MolFromSmarts("[#6]#[#7]")
    isoxazole_pattern = Chem.MolFromSmarts("[#6]1[#7][#8][#6][#7]1")

    # Track if we found the pattern
    found_pattern = False

    # Track the depth where nitrile disappears
    nitrile_disappears_at_depth = None

    # Track if isoxazole appears in final product
    isoxazole_in_final_product = False

    def dfs_traverse(node, depth=0):
        nonlocal found_pattern, nitrile_disappears_at_depth, isoxazole_in_final_product, findings_json

        if node["type"] == "mol":
            mol = Chem.MolFromSmiles(node["smiles"])
            if mol:
                # Check if this is the final product (depth 0)
                if depth == 0:
                    if mol.HasSubstructMatch(isoxazole_pattern):
                        isoxazole_in_final_product = True
                        if "isoxazole" not in findings_json["atomic_checks"]["ring_systems"]:
                            findings_json["atomic_checks"]["ring_systems"].append("isoxazole")

                # Check for nitrile presence
                has_nitrile = mol.HasSubstructMatch(nitrile_pattern)
                if has_nitrile:
                    if "nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("nitrile")

                # Process children (reactions)
                for child in node.get("children", []):
                    if child["type"] == "reaction":
                        # Get the product of this reaction (which is the current molecule)
                        product_smiles = node["smiles"]
                        product_mol = mol

                        # Get reactants
                        reactants = []
                        for reactant_node in child.get("children", []):
                            if reactant_node["type"] == "mol":
                                reactant_mol = Chem.MolFromSmiles(reactant_node["smiles"])
                                if reactant_mol:
                                    reactants.append(reactant_mol)

                        # Check if any reactant has nitrile but product doesn't
                        if not has_nitrile:
                            for reactant_mol in reactants:
                                if reactant_mol.HasSubstructMatch(nitrile_pattern):
                                    nitrile_disappears_at_depth = depth
                                    break

                        # Continue traversal
                        # Depth increases when going from chemical to reaction
                        dfs_traverse(child, depth + 1)
        else:  # node["type"] == "reaction"
            for child in node.get("children", []):
                # Depth remains the same when going from reaction to chemical
                dfs_traverse(child, depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we found the pattern: nitrile disappears in the final step (depth 0)
    # and isoxazole appears in the final product
    if nitrile_disappears_at_depth == 0 and isoxazole_in_final_product:
        found_pattern = True
        # Add structural constraints if the pattern is found
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "isoxazole",
                "position": "last_stage"
            }
        })
        findings_json["structural_constraints"].append({
            "type": "sequence",
            "details": {
                "before": "nitrile",
                "after": "isoxazole"
            }
        })
        if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
            findings_json["atomic_checks"]["named_reactions"].append("ring_formation")

    return found_pattern, findings_json
