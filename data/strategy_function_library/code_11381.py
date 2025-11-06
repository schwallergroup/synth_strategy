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
    This function detects a synthetic strategy involving protection of an amine
    as a carbamate, particularly a benzyl carbamate (Cbz).
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

    # Track if we found the protection
    found_protection = False

    def dfs_traverse(node, depth):
        nonlocal found_protection, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if not product or not all(reactants):
                print("Warning: Could not parse all molecules in reaction")
                return

            # Check for protection (amine to carbamate)
            carbamate_pattern = Chem.MolFromSmarts("[N][C](=[O])[O]")
            benzyl_carbamate_pattern = Chem.MolFromSmarts("[N][C](=[O])[O][CH2][c]")

            # Check if a carbamate is formed
            if (
                product.HasSubstructMatch(carbamate_pattern)
                and not any(r.HasSubstructMatch(carbamate_pattern) for r in reactants)
            ):
                found_protection = True
                findings_json["atomic_checks"]["named_reactions"].append("carbamate_formation")
                findings_json["atomic_checks"]["functional_groups"].append("carbamate")
                findings_json["structural_constraints"].append({
                    "type": "count",
                    "details": {
                        "target": "carbamate_formation",
                        "operator": ">=",
                        "value": 1
                    }
                })

                # Check if it's specifically a benzyl carbamate (Cbz)
                if product.HasSubstructMatch(benzyl_carbamate_pattern):
                    print("Found protection step: amine to benzyl carbamate (Cbz)")
                    findings_json["atomic_checks"]["functional_groups"].append("benzyl_carbamate")
                else:
                    print("Found protection step: amine to carbamate")

        # Process children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # Depth remains the same when traversing from reaction to chemical
                dfs_traverse(child, depth)
            else:
                # Depth increases when traversing from chemical to reaction
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route, 0) # Initial call with depth 0

    return found_protection, findings_json
