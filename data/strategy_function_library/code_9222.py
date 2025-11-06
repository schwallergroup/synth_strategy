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
    This function detects if the synthesis involves formation of a tricyclic core system.
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

    tricyclic_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal tricyclic_formation_detected, findings_json

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["mapped_reaction_smiles"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product_mol = Chem.MolFromSmiles(product_smiles)

            if product_mol and all(r for r in reactant_mols):
                # Count rings in reactants and product
                reactant_ring_counts = [len(Chem.GetSSSR(r)) for r in reactant_mols]
                product_ring_count = len(Chem.GetSSSR(product_mol))

                # Check if product has significantly more rings than any reactant
                if product_ring_count >= 3 and product_ring_count > max(reactant_ring_counts):
                    print(
                        f"Tricyclic formation detected: Product has {product_ring_count} rings, reactants have {reactant_ring_counts}"
                    )
                    tricyclic_formation_detected = True
                    # Record findings
                    if product_ring_count >= 3:
                        findings_json["structural_constraints"].append({
                            "type": "count",
                            "details": {
                                "target": "rings_in_product",
                                "operator": ">=",
                                "value": 3
                            }
                        })
                    # Assuming 'ring_formation' is a generic named reaction for this type of event
                    findings_json["atomic_checks"]["named_reactions"].append("ring_formation")
                    # Add the product ring system count as a finding
                    findings_json["atomic_checks"]["ring_systems"].append(f"{product_ring_count} rings")

        # Traverse children
        for child in node.get("children", []):
            if node["type"] == "reaction":
                dfs_traverse(child, depth)
            else:
                dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)
    return tricyclic_formation_detected, findings_json
