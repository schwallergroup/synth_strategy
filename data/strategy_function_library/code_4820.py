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
    Detects if a key ring-forming step occurs in the latter half of a synthesis (i.e., is 'late-stage').
    The function traverses the synthetic route to find the total number of steps (`max_depth`) and the step with the lowest depth value where a ring is formed.
    A ring formation is considered late-stage if its depth is less than or equal to half of the total synthesis depth (`depth=1` is the final synthetic step).
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

    # Track if we found a ring disconnection and at what depth
    found_disconnection = False
    disconnection_depth = float("inf")  # Initialize to infinity to find minimum depth
    max_depth = -1

    def dfs_traverse(node, depth=0):
        nonlocal found_disconnection, disconnection_depth, max_depth, findings_json

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["mapped_reaction_smiles"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                print(f"Analyzing reaction at depth {depth}:")
                print(f"  Product: {product_smiles}")
                print(f"  Reactants: {reactants_smiles}")

                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                if product and all(r for r in reactants):
                    # Count rings in product and reactants
                    product_rings = product.GetRingInfo().NumRings()
                    reactants_rings = sum(r.GetRingInfo().NumRings() for r in reactants)

                    print(f"  Ring count - Product: {product_rings}, Reactants: {reactants_rings}")

                    # Check for ring disconnection (product has more rings than reactants in retrosynthesis)
                    if product_rings > reactants_rings:
                        found_disconnection = True
                        findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                        # Update disconnection_depth if this is the latest stage (lowest depth)
                        if depth < disconnection_depth:
                            disconnection_depth = depth
                            print(f"  Found ring disconnection at depth {depth} (latest so far)")
            except Exception as e:
                print(f"Error processing reaction node: {e}")

        # Traverse children
        for child in node.get("children", []):
            # New logic for depth calculation
            if node["type"] == "reaction":
                # If current node is a reaction, depth remains the same for children
                dfs_traverse(child, depth)
            else:
                # If current node is not a reaction (e.g., chemical), depth increases
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # If no disconnection was found, set disconnection_depth to -1
    if not found_disconnection:
        disconnection_depth = -1
    else:
        # Add structural constraint for count if a ring formation was found
        findings_json["structural_constraints"].append({
            "type": "count",
            "details": {
                "target": "ring_formation",
                "operator": ">=",
                "value": 1
            }
        })

    # Check if the disconnection is in the second half of the synthesis
    # (low depth in retrosynthesis corresponds to late in synthesis)
    # For late stage disconnection, we want disconnection_depth to be low
    is_late_stage = found_disconnection and disconnection_depth <= max_depth / 2

    print(f"Disconnection depth: {disconnection_depth}, Max depth: {max_depth}")
    print(f"Is late stage: {is_late_stage}")

    if is_late_stage:
        # Add structural constraint for positional if it's late stage
        findings_json["structural_constraints"].append({
            "type": "positional",
            "details": {
                "target": "ring_formation",
                "position": "relative_late",
                "description": "The latest ring-forming step (one with the lowest depth value in the retrosynthesis tree) must occur in the latter half of the synthesis, defined as its depth being less than or equal to half the total route depth (depth <= max_depth / 2)."
            }
        })

    return is_late_stage, findings_json
