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
    This function detects if the synthesis route involves forming an ether linkage
    between two aromatic rings (specifically pyridine and nitrobenzene).
    """
    ether_formation_detected = False

    def dfs_traverse(node, depth=0):
        nonlocal ether_formation_detected

        if node["type"] == "reaction" and "rsmi" in node["metadata"]:
            rsmi = node["metadata"]["rsmi"]
            reactants = rsmi.split(">")[0].split(".")
            product = rsmi.split(">")[-1]

            reactant_mols = [
                Chem.MolFromSmiles(r) for r in reactants if Chem.MolFromSmiles(r)
            ]
            product_mol = Chem.MolFromSmiles(product)

            if product_mol and reactant_mols:
                # Check for ether linkage in product
                ether_pattern = Chem.MolFromSmarts("[c]-[#8]-[#6]")

                # Check for pyridine and nitrobenzene in product
                pyridine_pattern = Chem.MolFromSmarts("c1ccncc1")
                nitrobenzene_pattern = Chem.MolFromSmarts("c1ccc(cc1)[N+](=[O])[O-]")

                if (
                    product_mol.HasSubstructMatch(ether_pattern)
                    and product_mol.HasSubstructMatch(pyridine_pattern)
                    and product_mol.HasSubstructMatch(nitrobenzene_pattern)
                ):

                    # Check if reactants don't have the ether linkage
                    reactants_have_ether = any(
                        mol.HasSubstructMatch(ether_pattern)
                        for mol in reactant_mols
                        if mol
                    )

                    if not reactants_have_ether:
                        ether_formation_detected = True
                        print(
                            f"Detected ether linkage formation at depth {depth}: {rsmi}"
                        )

        for child in node.get("children", []):
            dfs_traverse(child, depth + 1)

    dfs_traverse(route)

    return ether_formation_detected
