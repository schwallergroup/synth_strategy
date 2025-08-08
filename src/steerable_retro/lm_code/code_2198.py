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
    Detects a synthesis strategy involving biaryl formation via Suzuki coupling.
    """
    # Track if we found Suzuki coupling for biaryl formation
    biaryl_suzuki_found = False

    def dfs_traverse(node, depth=0):
        nonlocal biaryl_suzuki_found

        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactants = [Chem.MolFromSmiles(r) for r in reactants_smiles if r]
            product = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if product and reactants:
                # Check for boronic acid/ester in reactants (indicator of Suzuki coupling)
                boronic_acid_pattern = Chem.MolFromSmarts("[#6]-[#5](-[#8])-[#8]")

                # Check for halide in reactants
                halide_pattern = Chem.MolFromSmarts("[#6]-[#9,#17,#35,#53]")  # C-F, C-Cl, C-Br, C-I

                # Check for biaryl in product
                # This is a simplified check - in reality would need more sophisticated analysis
                aromatic_pattern = Chem.MolFromSmarts("c:c")

                reactants_have_boronic = any(
                    r.HasSubstructMatch(boronic_acid_pattern) for r in reactants if r
                )
                reactants_have_halide = any(
                    r.HasSubstructMatch(halide_pattern) for r in reactants if r
                )
                product_has_aromatic = product.HasSubstructMatch(aromatic_pattern)

                # Simple heuristic for Suzuki coupling
                if reactants_have_boronic and reactants_have_halide and product_has_aromatic:
                    # Count aromatic rings in reactants and product
                    aromatic_rings_in_product = len(Chem.GetSSSR(product))
                    aromatic_rings_in_reactants = sum(len(Chem.GetSSSR(r)) for r in reactants if r)

                    # If product has same or more aromatic rings, might be biaryl formation
                    if aromatic_rings_in_product >= aromatic_rings_in_reactants:
                        biaryl_suzuki_found = True
                        print(f"Found potential biaryl Suzuki coupling at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    print(f"Biaryl formation via Suzuki coupling strategy detected: {biaryl_suzuki_found}")
    return biaryl_suzuki_found
