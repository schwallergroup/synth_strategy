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
    This function detects a strategy where stereochemistry is preserved
    throughout the synthesis (no changes to existing stereocenters).
    """
    preserves_stereochemistry = True

    def dfs_traverse(node, depth=0):
        nonlocal preserves_stereochemistry

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_part = rsmi.split(">")[0]
            product_part = rsmi.split(">")[-1]

            reactants = reactants_part.split(".")
            product = product_part

            # Convert to RDKit molecules
            try:
                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [
                    Chem.MolFromSmiles(r) for r in reactants if not r.startswith("CCOC")
                ]  # Ignore reagents

                if product_mol and len(reactant_mols) > 0:
                    main_reactant = reactant_mols[0]  # Assuming the first reactant is the main one

                    # Count chiral centers
                    reactant_chiral_centers = len(
                        Chem.FindMolChiralCenters(main_reactant, includeUnassigned=False)
                    )
                    product_chiral_centers = len(
                        Chem.FindMolChiralCenters(product_mol, includeUnassigned=False)
                    )

                    # If the number of chiral centers changes, stereochemistry is not preserved
                    if reactant_chiral_centers != product_chiral_centers:
                        print(
                            f"Stereochemistry changed at depth {depth}: {reactant_chiral_centers} → {product_chiral_centers} chiral centers"
                        )
                        preserves_stereochemistry = False
            except:
                print(f"Error processing reaction at depth {depth}")

        # Process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    if preserves_stereochemistry:
        print("Detected stereochemistry preservation strategy")

    return preserves_stereochemistry
