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
from typing import Tuple, Dict, List


def main(route) -> Tuple[bool, Dict]:
    """
    Detects if the synthetic route involves sequential modifications of the same functional group sites.
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

    phenol_modifications = 0
    lactam_n_modifications = 0

    def dfs_traverse(node, depth):
        nonlocal phenol_modifications, lactam_n_modifications, findings_json

        if node["type"] == "reaction" and "rsmi" in node.get("metadata", {}):
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_str = rsmi.split(">")[0]
            product_str = rsmi.split(">")[-1]

            reactants = [Chem.MolFromSmiles(r) for r in reactants_str.split(".") if r]
            product = Chem.MolFromSmiles(product_str)

            if product and all(r for r in reactants):
                # Check for lactam N modifications
                lactam_n_pattern = Chem.MolFromSmarts("[NH][C]=[O]")
                n_modified_pattern = Chem.MolFromSmarts("[N]([#6])[C]=[O]")

                # Check if a reactant has a lactam N and the product has a modified lactam N
                if any(
                    r.HasSubstructMatch(lactam_n_pattern) for r in reactants if r
                ) and product.HasSubstructMatch(n_modified_pattern):
                    lactam_n_modifications += 1
                    # Record atomic checks for lactam N modification
                    if "secondary_amide_or_lactam" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("secondary_amide_or_lactam")
                    if "tertiary_amide_or_lactam" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("tertiary_amide_or_lactam")
                    if "lactam_N_modification" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("lactam_N_modification")

                # Check for phenol modifications (conceptual, as original code didn't have explicit phenol check)
                # For demonstration, let's assume a placeholder check for phenol modification
                # In a real scenario, this would involve specific SMARTS patterns or reaction types
                # For now, we'll just increment if a certain condition (e.g., a specific reaction type) is met
                # As the original code only increments phenol_modifications without a specific check,
                # we'll assume it's handled elsewhere or is a placeholder for future implementation.
                # To make it functional for this refactoring, let's add a dummy condition.
                # If the original code had a specific check, it would be placed here.
                # For the purpose of this exercise, we'll assume a condition that would increment phenol_modifications.
                # Since the original code doesn't have a specific check for phenol modification within dfs_traverse,
                # we cannot add an atomic check for it here based on the provided code.
                # The `phenol_modifications` counter is present, but its increment logic is missing in the provided snippet.
                # We will only record findings for `lactam_n_modifications` as its logic is present.

        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    dfs_traverse(route, 0)

    # Strategy is present if there are multiple modifications of at least one functional group
    result = phenol_modifications >= 2 or lactam_n_modifications >= 2

    # Record structural constraints if met
    if phenol_modifications >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "phenol_modification",
                "operator": ">=",
                "value": 2
            }
        })
    if lactam_n_modifications >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "lactam_N_modification",
                "operator": ">=",
                "value": 2
            }
        })

    print(
        f"Sequential functional group modifications detected: {result} (phenol: {phenol_modifications}, lactam N: {lactam_n_modifications})"
    )
    return result, findings_json