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
    Detects a synthetic strategy involving nitro reduction followed by amide coupling
    in a linear synthesis pathway.
    """
    # Track transformations and their depths
    nitro_reduction_depth = None
    amide_formation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal nitro_reduction_depth, amide_formation_depth

        if node["type"] == "reaction":
            # Extract reactants and product from reaction SMILES
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if r]
            product_mol = Chem.MolFromSmiles(product) if product else None

            if not product_mol or not all(reactant_mols):
                print("Warning: Could not parse some molecules in reaction")
                return

            # Check for nitro reduction
            nitro_pattern = Chem.MolFromSmarts("[N+](=[O])[O-]")
            amine_pattern = Chem.MolFromSmarts("[NH2]")

            if any(
                mol.HasSubstructMatch(nitro_pattern) for mol in reactant_mols
            ) and product_mol.HasSubstructMatch(amine_pattern):
                nitro_reduction_depth = depth
                print(f"Found nitro reduction at depth {depth}")

            # Check for amide formation
            amine_reactant = any(mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols)
            amide_pattern = Chem.MolFromSmarts("[NH]C(=O)")
            acyl_halide_pattern = Chem.MolFromSmarts("C(=O)[Cl,Br,I]")

            acyl_halide_present = any(
                mol.HasSubstructMatch(acyl_halide_pattern) for mol in reactant_mols
            )

            if (
                amine_reactant
                and acyl_halide_present
                and product_mol.HasSubstructMatch(amide_pattern)
            ):
                amide_formation_depth = depth
                print(f"Found amide formation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal from the root
    dfs_traverse(route)

    # Check if we found both transformations and if nitro reduction comes before amide formation
    strategy_present = (
        nitro_reduction_depth is not None
        and amide_formation_depth is not None
        and nitro_reduction_depth > amide_formation_depth
    )  # Higher depth means earlier in synthesis

    print(f"Strategy detection result: {strategy_present}")
    return strategy_present
