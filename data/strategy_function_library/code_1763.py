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
    This function detects if the synthetic route involves nitrile-mediated heterocycle formation,
    where a nitrile group is maintained through multiple steps before being used in cyclization.
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

    nitrile_present_steps = 0
    nitrile_used_in_cyclization = False
    result = False

    def dfs_traverse(node, depth=0):
        nonlocal nitrile_present_steps, nitrile_used_in_cyclization, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols if r):
                # Check for nitrile presence
                nitrile_pattern = Chem.MolFromSmarts("[C]#[N]")
                pyrazole_pattern = Chem.MolFromSmarts("c1[n]n[c]c1")

                reactants_have_nitrile = any(
                    r.HasSubstructMatch(nitrile_pattern) for r in reactant_mols if r
                )
                product_has_nitrile = product_mol.HasSubstructMatch(nitrile_pattern)

                # Count steps where nitrile is present
                if reactants_have_nitrile or product_has_nitrile:
                    nitrile_present_steps += 1
                    if "nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                        findings_json["atomic_checks"]["functional_groups"].append("nitrile")

                # Check if nitrile is used in heterocycle formation
                if (
                    reactants_have_nitrile
                    and not product_has_nitrile
                    and product_mol.HasSubstructMatch(pyrazole_pattern)
                ):
                    nitrile_used_in_cyclization = True
                    if "pyrazole" not in findings_json["atomic_checks"]["ring_systems"]:
                        findings_json["atomic_checks"]["ring_systems"].append("pyrazole")
                    # This implies a ring formation reaction
                    if "ring_formation" not in findings_json["atomic_checks"]["named_reactions"]:
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    findings_json["structural_constraints"].append({
                        "type": "co-occurrence",
                        "details": {
                            "scope": "single_reaction",
                            "targets": [
                                "nitrile_in_reactant",
                                "pyrazole_in_product",
                                "negation(nitrile_in_product)"
                            ]
                        }
                    })

        # Traverse children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                # If current node is reaction, depth remains the same for child (chemical)
                dfs_traverse(child, depth)
            else:
                # If current node is chemical, depth increases for child (reaction)
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine final result and update findings_json for count constraint
    result = nitrile_present_steps >= 2 and nitrile_used_in_cyclization

    if nitrile_present_steps >= 2:
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "reaction_with_nitrile",
                "operator": ">=",
                "value": 2
            }
        })

    # Return True if nitrile is present in multiple steps and used in cyclization
    return result, findings_json
