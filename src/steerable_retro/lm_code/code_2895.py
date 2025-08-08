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
    This function detects a linear synthesis strategy with sequential aromatic functionalization.
    """
    # Track key transformations in sequence
    transformations = []

    def dfs_traverse(node, depth=0):
        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0].split(".")
                product = rsmi.split(">")[-1]

                # Analyze the type of transformation
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]
                product_mol = Chem.MolFromSmiles(product)

                # Check for aromatic core
                aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")
                has_aromatic = product_mol and product_mol.HasSubstructMatch(aromatic_pattern)

                if has_aromatic:
                    # Determine transformation type
                    if any(
                        mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[c]-[Br,Cl,I]"))
                        for mol in reactant_mols
                    ):
                        transformations.append(("aryl_coupling", depth))
                    elif any(
                        mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[N+](=[O])[O-]"))
                        for mol in reactant_mols
                    ):
                        transformations.append(("nitro_reduction", depth))
                    elif any(
                        mol and mol.HasSubstructMatch(Chem.MolFromSmarts("[C](=[O])[O]"))
                        for mol in reactant_mols
                    ):
                        transformations.append(("amide_formation", depth))

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    # Check if we have a linear sequence of aromatic functionalizations
    if len(transformations) >= 2:
        # Sort by depth (higher depth = earlier in synthesis)
        transformations.sort(key=lambda x: x[1], reverse=True)
        transformation_types = [t[0] for t in transformations]

        print(f"Detected transformation sequence: {transformation_types}")

        # Check for specific patterns that indicate linear aromatic functionalization
        if "aryl_coupling" in transformation_types and (
            "nitro_reduction" in transformation_types or "amide_formation" in transformation_types
        ):
            return True

    return False
