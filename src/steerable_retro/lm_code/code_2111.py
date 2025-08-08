#!/bin/python

"""LM-defined function for strategy description."""

from rdkit.Chem import AllChem, rdFMCS
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


def main(route):
    """
    This function detects a linear synthesis strategy with multiple ring transformations
    (at least one ring formation and one ring opening).
    """
    # Track ring transformations
    ring_formations = 0
    ring_openings = 0

    def dfs_traverse(node):
        nonlocal ring_formations, ring_openings

        if node["type"] == "reaction":
            try:
                # Extract reactants and product
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0].split(".")
                product_smiles = rsmi.split(">")[-1]

                # Convert to RDKit molecules
                reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product = Chem.MolFromSmiles(product_smiles)

                # Check if all molecules were parsed correctly
                if product and all(r is not None for r in reactants):
                    # Count rings in reactants and product
                    product_ring_count = product.GetRingInfo().NumRings()
                    reactant_ring_count = sum(r.GetRingInfo().NumRings() for r in reactants)

                    if product_ring_count > reactant_ring_count:
                        print(f"Detected ring formation: {rsmi}")
                        print(
                            f"  Product rings: {product_ring_count}, Reactant rings: {reactant_ring_count}"
                        )
                        ring_formations += 1
                    elif product_ring_count < reactant_ring_count:
                        print(f"Detected ring opening: {rsmi}")
                        print(
                            f"  Product rings: {product_ring_count}, Reactant rings: {reactant_ring_count}"
                        )
                        ring_openings += 1
            except Exception as e:
                print(f"Error processing reaction: {e}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if there's at least one ring formation and one ring opening
    result = ring_formations >= 1 and ring_openings >= 1
    print(f"Ring formations: {ring_formations}, Ring openings: {ring_openings}")
    return result
