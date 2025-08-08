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
    This function detects a strategy where a sulfonamide group is involved in both
    ring opening and subsequent ring formation steps.
    """
    # Track reactions and their properties
    reactions_with_sulfonamide = []

    # SMARTS patterns
    sulfonamide_pattern = "[#7][#16](=[#8])(=[#8])[#6]"

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check if sulfonamide is involved
            if has_substructure(product, sulfonamide_pattern) or any(
                has_substructure(r, sulfonamide_pattern) for r in reactants
            ):

                # Count rings before and after
                reactant_rings = sum(count_rings(r) for r in reactants)
                product_rings = count_rings(product)
                ring_change = product_rings - reactant_rings

                reactions_with_sulfonamide.append(
                    {"depth": depth, "ring_change": ring_change, "rsmi": rsmi}
                )
                print(f"Found sulfonamide reaction at depth {depth} with ring change {ring_change}")

        # Recursively process children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    def has_substructure(smiles, pattern):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            patt = Chem.MolFromSmarts(pattern)
            return mol.HasSubstructMatch(patt)
        return False

    def count_rings(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return len(Chem.GetSSSR(mol))
        return 0

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we have both ring opening and formation with sulfonamide involvement
    has_ring_opening = any(r["ring_change"] < 0 for r in reactions_with_sulfonamide)
    has_ring_formation = any(r["ring_change"] > 0 for r in reactions_with_sulfonamide)

    return has_ring_opening and has_ring_formation
