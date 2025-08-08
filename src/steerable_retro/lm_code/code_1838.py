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
    This function detects a convergent synthesis with late-stage amide formation.
    It looks for amide bond formation in the second half of the synthesis
    and checks if multiple fragments are combined.
    """
    # Track if we found amide formation and at what depth
    amide_formation_found = False
    amide_formation_depth = -1
    fragment_count = 0
    max_depth = 0

    def dfs_traverse(node, depth=0):
        nonlocal amide_formation_found, amide_formation_depth, fragment_count, max_depth

        max_depth = max(max_depth, depth)

        if node["type"] == "reaction":
            # Check if this is an amide formation reaction
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Check for acid chloride in reactants
                acid_chloride_present = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[Cl]"))
                    for r in reactants
                    if r
                )

                # Check for amine in reactants
                amine_present = any(
                    Chem.MolFromSmiles(r) is not None
                    and Chem.MolFromSmiles(r).HasSubstructMatch(Chem.MolFromSmarts("[N;H2]"))
                    for r in reactants
                    if r
                )

                # Check for amide in product
                product_mol = Chem.MolFromSmiles(product) if product else None
                amide_in_product = product_mol is not None and product_mol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C](=[O])[N]")
                )

                # If we have acid chloride + amine → amide, it's an amide formation
                if acid_chloride_present and amine_present and amide_in_product:
                    amide_formation_found = True
                    amide_formation_depth = depth
                    print(f"Found amide formation at depth {depth}")

                # Count fragments being combined
                if len(reactants) > 1:
                    fragment_count += len(reactants) - 1
                    print(f"Found {len(reactants)} fragments being combined at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Determine if this is a late-stage amide formation in a convergent synthesis
    is_late_stage = amide_formation_depth >= 0 and amide_formation_depth <= max_depth / 2
    is_convergent = fragment_count >= 2

    print(
        f"Amide formation: {amide_formation_found}, Depth: {amide_formation_depth}, Max depth: {max_depth}"
    )
    print(f"Fragment count: {fragment_count}")
    print(f"Is late stage: {is_late_stage}, Is convergent: {is_convergent}")

    return amide_formation_found and is_late_stage and is_convergent
