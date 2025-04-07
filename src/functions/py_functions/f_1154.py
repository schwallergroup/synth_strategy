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
    Detects if the synthesis involves multiple SNAr reactions for C-N bond formation,
    particularly on pyrimidine or similar electron-deficient aromatic rings.
    """
    snar_count = 0

    def dfs_traverse(node):
        nonlocal snar_count

        if (
            node.get("type") == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for SNAr pattern: aromatic carbon with leaving group (Cl, F, Br)
            # being replaced by nitrogen nucleophile

            # Create molecules
            product_mol = Chem.MolFromSmiles(product)
            reactant_mols = [Chem.MolFromSmiles(r) for r in reactants]

            if product_mol and all(reactant_mols):
                # Look for pyrimidine or similar electron-deficient aromatic rings in reactants
                pyrimidine_pattern = Chem.MolFromSmarts("c1ncncc1")

                # Check if any reactant contains pyrimidine
                has_pyrimidine = False
                for r_mol in reactant_mols:
                    if r_mol.HasSubstructMatch(pyrimidine_pattern):
                        has_pyrimidine = True
                        break

                # Look for C-N bond formation
                # This is a simplified check - in a real implementation, you would need to
                # analyze the reaction more carefully using atom mapping
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                halide_pattern = Chem.MolFromSmarts("c[Cl,F,Br]")

                has_amine = False
                has_halide = False

                for r_mol in reactant_mols:
                    if r_mol.HasSubstructMatch(amine_pattern):
                        has_amine = True
                    if r_mol.HasSubstructMatch(halide_pattern):
                        has_halide = True

                # If we have a pyrimidine, an amine, and a halide, it's likely an SNAr
                if has_pyrimidine and has_amine and has_halide:
                    snar_count += 1
                    print(f"SNAr reaction detected: {rsmi}")

        # Traverse children
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal
    dfs_traverse(route)

    # Check if we have multiple SNAr reactions
    has_multiple_snar = snar_count >= 2
    print(f"Number of SNAr reactions detected: {snar_count}")

    return has_multiple_snar
