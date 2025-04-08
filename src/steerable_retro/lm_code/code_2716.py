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
    Detects if the synthesis involves the formation of a tetracyclic system
    through a cyclization reaction.
    """
    tetracyclic_formation_detected = False

    def dfs_traverse(node):
        nonlocal tetracyclic_formation_detected

        if node["type"] == "reaction":
            if "rsmi" in node.get("metadata", {}):
                rsmi = node["metadata"]["rsmi"]
                reactants = rsmi.split(">")[0]
                product = rsmi.split(">")[-1]

                product_mol = Chem.MolFromSmiles(product)
                reactant_mols = [Chem.MolFromSmiles(r) for r in reactants.split(".")]

                if product_mol and all(reactant_mol for reactant_mol in reactant_mols):
                    # Count rings in product and reactants
                    product_rings = Chem.GetSSSR(product_mol)

                    max_reactant_rings = 0
                    for r_mol in reactant_mols:
                        if r_mol:
                            r_rings = Chem.GetSSSR(r_mol)
                            max_reactant_rings = max(max_reactant_rings, len(r_rings))

                    # Check if product has at least 4 rings and more rings than any reactant
                    if len(product_rings) >= 4 and len(product_rings) > max_reactant_rings:
                        tetracyclic_formation_detected = True
                        print(
                            f"Tetracyclic formation detected: Product has {len(product_rings)} rings, max reactant has {max_reactant_rings} rings"
                        )

        for child in node.get("children", []):
            dfs_traverse(child)

    dfs_traverse(route)
    return tetracyclic_formation_detected
