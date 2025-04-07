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
    Detects if the synthesis route involves a Suzuki coupling reaction,
    looking for boronic acid/ester reactants and C-C bond formation.
    """
    result = False

    def dfs_traverse(node):
        nonlocal result

        if (
            node["type"] == "reaction"
            and "metadata" in node
            and "rsmi" in node["metadata"]
        ):
            rsmi = node["metadata"]["rsmi"]
            reactants_smiles = rsmi.split(">")[0].split(".")
            product_smiles = rsmi.split(">")[-1]

            try:
                # Parse molecules
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants_smiles]
                product_mol = Chem.MolFromSmiles(product_smiles)

                if product_mol and all(r for r in reactant_mols):
                    # Check for boronic acid/ester pattern in reactants
                    boronic_pattern = Chem.MolFromSmarts("[c,C]-[B]([O,OH])[O,OH]")

                    # Check for aryl halide pattern in reactants
                    aryl_halide_pattern = Chem.MolFromSmarts("c-[Cl,Br,I]")

                    if any(
                        r.HasSubstructMatch(boronic_pattern) for r in reactant_mols
                    ) and any(
                        r.HasSubstructMatch(aryl_halide_pattern) for r in reactant_mols
                    ):
                        # This is a potential Suzuki coupling
                        print("Detected Suzuki coupling reaction")
                        result = True

            except Exception as e:
                print(f"Error in suzuki_coupling_in_synthesis: {e}")

        # Continue traversal
        for child in node.get("children", []):
            dfs_traverse(child)

    # Start traversal from the root
    dfs_traverse(route)
    return result
