#!/bin/python

"""LM-defined function for strategy description."""

import copy
import re
from collections import deque

import rdkit
import rdkit.Chem as Chem
from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Lipinski,
    rdChemReactions,
    rdFMCS,
    rdMolDescriptors,
    rdmolops,
)
from rdkit.Chem.Scaffolds import MurckoScaffold


def main(route):
    """
    This function detects a strategy involving S-alkylation following sulfonylation.
    """
    # Track reactions and their depths
    sulfonylation_depth = None
    s_alkylation_depth = None

    def dfs_traverse(node, depth=0):
        nonlocal sulfonylation_depth, s_alkylation_depth

        # Only process reaction nodes
        if node["type"] == "reaction":
            # Extract reactants and product
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            # Convert to RDKit molecules
            reactant_mols = [Chem.MolFromSmiles(smi) for smi in reactants_smiles if smi]
            product_mol = Chem.MolFromSmiles(product_smiles) if product_smiles else None

            if all(reactant_mols) and product_mol:
                # Check for sulfonylation (formation of sulfonyl chloride)
                sulfonyl_cl_pattern = Chem.MolFromSmarts("[c:1][S:2](=[O:3])(=[O:4])[Cl:5]")
                if product_mol.HasSubstructMatch(sulfonyl_cl_pattern):
                    has_sulfonyl_cl_reactant = any(
                        mol.HasSubstructMatch(sulfonyl_cl_pattern) for mol in reactant_mols if mol
                    )
                    if not has_sulfonyl_cl_reactant:
                        sulfonylation_depth = depth
                        print(f"Found sulfonylation at depth {depth}")

                # Check for S-alkylation (conversion of SO2Cl to SO2R)
                s_alkyl_pattern = Chem.MolFromSmarts("[c:1][S:2](=[O:3])(=[O:4])[#6:5]")
                sulfonyl_cl_pattern = Chem.MolFromSmarts("[c:1][S:2](=[O:3])(=[O:4])[Cl:5]")

                if product_mol.HasSubstructMatch(s_alkyl_pattern):
                    # Check if any reactant has sulfonyl chloride
                    has_sulfonyl_cl_reactant = any(
                        mol and mol.HasSubstructMatch(sulfonyl_cl_pattern) for mol in reactant_mols
                    )
                    if has_sulfonyl_cl_reactant:
                        s_alkylation_depth = depth
                        print(f"Found S-alkylation at depth {depth}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    # Start traversal
    dfs_traverse(route)

    # Return True if S-alkylation follows sulfonylation
    if sulfonylation_depth is not None and s_alkylation_depth is not None:
        return (
            s_alkylation_depth < sulfonylation_depth
        )  # Remember: lower depth = later in synthesis
    return False
