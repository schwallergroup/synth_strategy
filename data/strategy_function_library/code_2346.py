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
    This function detects the use of Weinreb amide for controlled reduction to aldehyde
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

    # Define SMARTS patterns
    weinreb_amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[N][OX2][#6]")
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=[OX1])[#6]")

    # Track if we find the patterns and at what depth
    found_weinreb_depth = None
    found_aldehyde_depth = None

    result = False

    def dfs_traverse(node, depth=0):
        nonlocal found_weinreb_depth, found_aldehyde_depth, findings_json

        # Check for Weinreb amide to aldehyde transformation in reaction nodes
        if node["type"] == "reaction" and "metadata" in node and "mapped_reaction_smiles" in node["metadata"]:
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mol = Chem.MolFromSmiles(reactants[0]) if reactants else None
            product_mol = Chem.MolFromSmiles(product) if product else None

            if reactant_mol and product_mol:
                if reactant_mol.HasSubstructMatch(
                    weinreb_amide_pattern
                ) and product_mol.HasSubstructMatch(aldehyde_pattern):
                    found_weinreb_depth = depth
                    found_aldehyde_depth = depth - 1  # Aldehyde appears one step later in synthesis
                    findings_json["atomic_checks"]["functional_groups"].append("Weinreb amide")
                    findings_json["atomic_checks"]["functional_groups"].append("aldehyde")
                    findings_json["atomic_checks"]["named_reactions"].append("Weinreb aldehyde synthesis")

        for child in node.get("children", []):
            # New logic for depth calculation
            new_depth = depth
            if node["type"] != "reaction": # If current node is chemical, depth increases
                new_depth = depth + 1
            # If current node is reaction, depth remains the same
            dfs_traverse(child, new_depth)

    # Start traversal
    dfs_traverse(route)

    # Check if we found both Weinreb amide and aldehyde in the correct sequence
    if found_weinreb_depth is not None and found_aldehyde_depth is not None:
        if found_weinreb_depth > found_aldehyde_depth:  # Weinreb comes before aldehyde in synthesis
            print(
                f"Found Weinreb amide strategy: Weinreb at depth {found_weinreb_depth}, aldehyde at depth {found_aldehyde_depth}"
            )
            result = True
            # Add the structural constraint if the main condition is met
            findings_json["structural_constraints"].append({
                "type": "co-occurrence",
                "details": {
                    "targets": [
                        "Weinreb aldehyde synthesis"
                    ]
                }
            })

    if not result:
        print("Did not find evidence of Weinreb amide strategy")

    return result, findings_json