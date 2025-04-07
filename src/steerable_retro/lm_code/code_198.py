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
    Detects if the synthesis involves fragment coupling via amide bond formation.
    """
    amide_formation_found = False

    def dfs_traverse(node):
        nonlocal amide_formation_found

        if node["type"] == "reaction" and "metadata" in node and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            # Check for amide formation pattern
            if len(reactants) >= 2:  # At least two reactants for fragment coupling
                # Look for acid chloride and amine patterns in reactants
                acid_chloride_pattern = Chem.MolFromSmarts("C(=O)Cl")
                amine_pattern = Chem.MolFromSmarts("[NH2]")
                amide_pattern = Chem.MolFromSmarts("C(=O)N")

                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)]
                product_mol = Chem.MolFromSmiles(product)

                if product_mol and product_mol.HasSubstructMatch(amide_pattern):
                    has_acid_chloride = any(
                        mol and mol.HasSubstructMatch(acid_chloride_pattern)
                        for mol in reactant_mols
                    )
                    has_amine = any(
                        mol and mol.HasSubstructMatch(amine_pattern) for mol in reactant_mols
                    )

                    if has_acid_chloride and has_amine:
                        amide_formation_found = True
                        print(
                            f"Found amide bond formation at depth {node.get('metadata', {}).get('ID', '')}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return amide_formation_found
