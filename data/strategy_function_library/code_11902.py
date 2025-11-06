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
    Detects if the synthesis route includes a late-stage amidoxime to nitrile conversion.
    Late stage means in the first half of the synthesis (low depth in retrosynthetic tree).
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

    found_conversion = False
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal found_conversion, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction" and depth <= 1:  # Late stage (first or second reaction)
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                # Check for amidoxime pattern in reactants and nitrile in product
                amidoxime_pattern = Chem.MolFromSmarts("[CX3](=[NX2])[NX3][OX2]")
                nitrile_pattern = Chem.MolFromSmarts("[CX2]#[NX1]")

                product_mol = Chem.MolFromSmiles(product)
                reactants_mol = Chem.MolFromSmiles(reactants)

                if product_mol and reactants_mol:
                    amidoxime_found_in_reactants = reactants_mol.HasSubstructMatch(amidoxime_pattern)
                    nitrile_found_in_product = product_mol.HasSubstructMatch(nitrile_pattern)

                    if amidoxime_found_in_reactants:
                        if "amidoxime" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("amidoxime")
                    if nitrile_found_in_product:
                        if "nitrile" not in findings_json["atomic_checks"]["functional_groups"]:
                            findings_json["atomic_checks"]["functional_groups"].append("nitrile")

                    if amidoxime_found_in_reactants and nitrile_found_in_product:
                        found_conversion = True
                        if "amidoxime_to_nitrile_conversion" not in findings_json["atomic_checks"]["named_reactions"]:
                            findings_json["atomic_checks"]["named_reactions"].append("amidoxime_to_nitrile_conversion")

        for child in node.get("children", []):
            # New logic for depth calculation:
            # Depth increases only when traversing from a chemical node to a reaction node.
            # Depth remains the same when traversing from a reaction node to a chemical node.
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:  # Assuming 'chemical' type or any other type that should increase depth
                dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    result = False
    # Check if conversion was found in the first half of synthesis
    if found_conversion and max_depth > 0:
        result = True
        # Add the structural constraint if the condition is met
        findings_json["structural_constraints"].append(
            {
                "type": "positional",
                "details": {
                    "target": "amidoxime_to_nitrile_conversion",
                    "position": "last_two_stages"
                }
            }
        )

    return result, findings_json
