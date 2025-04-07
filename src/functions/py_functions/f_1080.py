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
    This function detects if the synthetic route involves manipulation of halogen atoms,
    specifically looking for introduction or removal of iodine while maintaining fluorine atoms.
    """
    aryl_iodide_pattern = Chem.MolFromSmarts("c[I]")
    has_halogen_manipulation = False

    def dfs_traverse(node):
        nonlocal has_halogen_manipulation

        if node["type"] == "reaction":
            if "metadata" in node and "rsmi" in node["metadata"]:
                rsmi = node["metadata"]["rsmi"]
                reactants_smiles = rsmi.split(">")[0]
                product_smiles = rsmi.split(">")[-1]

                try:
                    reactants_mol = Chem.MolFromSmiles(reactants_smiles)
                    product_mol = Chem.MolFromSmiles(product_smiles)

                    # Check for iodine addition or removal
                    if reactants_mol and product_mol:
                        reactant_has_iodine = reactants_mol.HasSubstructMatch(
                            aryl_iodide_pattern
                        )
                        product_has_iodine = product_mol.HasSubstructMatch(
                            aryl_iodide_pattern
                        )

                        if reactant_has_iodine != product_has_iodine:
                            print(f"Detected halogen manipulation in reaction: {rsmi}")
                            has_halogen_manipulation = True
                except:
                    print(f"Error processing reaction SMILES: {rsmi}")

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return has_halogen_manipulation
